process JUICER {

    tag "$sample_id"

    publishDir "${params.outdir}/hic", mode: 'copy'

    input:
    tuple val(sample_id), path(pairs)

    output:
    path "${sample_id}.hic"

    script:
    """
    juicer_tools pre \
        ${pairs} \
        ${sample_id}.hic \
        ${params.chromsz}
    """
}
