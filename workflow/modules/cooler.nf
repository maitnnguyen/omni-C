process COOLER {
    tag "${sample_id}"
    label 'process_high'

    publishDir "${params.outdir}/cooler/${sample_id}", mode: 'copy'

    input:
    tuple val(sample_id), path(pairs)
    path  chrom_sizes

    output:
    tuple val(sample_id), path("${sample_id}.mcool"), emit: mcool

    script:
    def prefix     = "${sample_id}"
    // Base resolution = finest resolution; must not be repeated in zoomify list.
    // Pull from params so users can override with --base_resolution 500 etc.
    def base_res   = params.base_resolution  ?: 1000
    // Zoomify resolutions — all coarser than base_res.
    // Kept in params so users can add/remove without touching module code.
    def zoom_res   = params.zoom_resolutions ?: "5000,10000,25000,50000,100000,500000,1000000"
    """
    # 1. Load pairs into a single-resolution .cool at the base resolution.
    #    --assembly informs cooler of the genome build (stored in .cool metadata).
    #    Cooler reads column positions from the pairs header automatically.
    cooler cload pairs \\
        --assembly ${params.genome} \\
        --chunksize 10000000 \\
        ${chrom_sizes}:${base_res} \\
        ${pairs} \\
        ${prefix}.base.cool

    # 2. Zoomify: build the multi-resolution .mcool from the base cool.
    #    --balance applies ICE balancing at every resolution during zoomify —
    #    this is faster than running cooler balance separately because cooler
    #    can reuse the aggregated matrix already in memory at each zoom level.
    #    --balance-args passes quality cutoffs to avoid unstable sparse matrices:
    #      --min-nnz 10   : skip rows/cols with fewer than 10 non-zero contacts
    #      --mad-max 5    : drop outlier bins more than 5 MADs from median
    cooler zoomify \\
        --nproc       ${task.cpus} \\
        --resolutions ${zoom_res} \\
        --balance \\
        --balance-args "--min-nnz 10 --mad-max 5 --nproc ${task.cpus}" \\
        -o ${prefix}.mcool \\
        ${prefix}.base.cool

    # Clean up base cool — it's fully embedded in the .mcool already
    rm ${prefix}.base.cool
    """

    stub:
    """
    touch ${sample_id}.mcool
    """
}