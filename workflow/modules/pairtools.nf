process PAIRTOOLS {

    tag "$sample_id"

    publishDir "${params.outdir}/pairs", mode: 'copy'

    input:
    tuple val(sample_id), path(bam)

    output:
    tuple val(sample_id), path("${sample_id}.pairs")

    script:
    """
    pairtools parse \
        --chroms-path ${params.chromsz} \
        ${bam} > ${sample_id}.pairs
    """
}
