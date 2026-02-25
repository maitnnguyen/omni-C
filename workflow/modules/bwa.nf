process ALIGN {

    tag "$sample_id"

    publishDir "${params.outdir}/alignment", mode: 'copy'

    input:
    tuple val(sample_id), path(r1), path(r2)

    output:
    tuple val(sample_id), path("${sample_id}.bam")

    script:
    """
    bwa mem ${params.genome} ${r1} ${r2} |
        samtools sort -o ${sample_id}.bam
    """
}
