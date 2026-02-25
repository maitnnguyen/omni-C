process COOLER {

    tag "$sample_id"

    publishDir "${params.outdir}/cooler", mode: 'copy'

    input:
    tuple val(sample_id), path(pairs)

    output:
    path "${sample_id}.cool"

    script:
    """
    cooler cload pairs \
        ${params.chromsz}:${params.binsize} \
        ${pairs} \
        ${sample_id}.cool
    """
}
