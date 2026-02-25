process FASTQC {

    tag "$sample_id"

    publishDir "${params.outdir}/fastqc", mode: 'copy'

    input:
    tuple val(sample_id), path(r1), path(r2)

    output:
    path "*.zip"
    path "*.html"

    script:
    """
    fastqc ${r1} ${r2}
    """
}
