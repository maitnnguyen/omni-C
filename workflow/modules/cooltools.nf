// modules/cooltools.nf

process COOLTOOLS_SCALING {
    tag "${sample_id} (${params.genome})"
    label 'process_medium'

    publishDir "${params.outdir}/qc/${sample_id}", mode: 'copy'

    input:
    tuple val(sample_id), path(mcool)
    path  chrom_sizes

    output:
    tuple val(sample_id), path("${sample_id}.expected.tsv"), emit: expected
    path  "${sample_id}.ps_curve.tsv",                        emit: ps_curve
    path  "${sample_id}.ps_curve.png",                        emit: plot

    script:
    def prefix  = "${sample_id}"
    def qc_res  = params.qc_resolution ?: 10000
    def uri     = "${mcool}::resolutions/${qc_res}"
    """
    # 1. Compute expected cis contacts (genome-wide, per-chromosome)
    cooltools expected-cis \\
        --nproc ${task.cpus} \\
        --clr-weight-name weight \\
        -o ${prefix}.expected.tsv \\
        ${uri}

    # 2. Compute P(s) curve from the cooler directly
    cooltools ps \\
        --nproc ${task.cpus} \\
        -o ${prefix}.ps_curve.tsv \\
        ${uri} \\
        ${uri}

    # 3. Plot P(s) in a dedicated Python script — not inline.
    #    The script is written here so the module is self-contained,
    #    but it could equally live in bin/ and be called by name.
    cat > plot_ps.py << 'EOF'
        import sys
        import pandas as pd
        import matplotlib
        matplotlib.use('Agg')          # Non-interactive backend — no display needed
        import matplotlib.pyplot as plt
        import numpy as np

        sample_id = sys.argv[1]
        ps_file   = sys.argv[2]
        out_png   = sys.argv[3]
        qc_res    = int(sys.argv[4])

        df = pd.read_csv(ps_file, sep='\\t', comment='#')

        # Drop empty bins and very short distances (sub-resolution noise)
        df = df[(df['balanced.avg'] > 0) & (df['dist'] >= qc_res)]

        fig, ax = plt.subplots(figsize=(6, 5))

        ax.loglog(df['dist'], df['balanced.avg'],
                color='steelblue', linewidth=1.5, label=sample_id)

        # Reference slope: expected power-law decay for a well-formed OmniC library
        # is approximately -1.5 in the 10kb–10Mb range
        x_ref = np.array([1e4, 1e7])
        y_ref = 1e-3 * (x_ref / 1e4) ** -1.5
        ax.loglog(x_ref, y_ref,
                color='gray', linewidth=1, linestyle='--', label='slope = -1.5')

        ax.set_xlabel('Genomic separation (bp)')
        ax.set_ylabel('Contact probability P(s)')
        ax.set_title(f'P(s) curve — {sample_id}')
        ax.legend(frameon=False)
        ax.grid(True, which='both', alpha=0.3)

        plt.tight_layout()
        plt.savefig(out_png, dpi=150)
        plt.close()
    EOF

    python plot_ps.py \\
        "${prefix}" \\
        "${prefix}.ps_curve.tsv" \\
        "${prefix}.ps_curve.png" \\
        "${params.qc_resolution ?: 10000}"
    """

    stub:
    """
    touch ${sample_id}.expected.tsv
    touch ${sample_id}.ps_curve.tsv
    touch ${sample_id}.ps_curve.png
    """
}