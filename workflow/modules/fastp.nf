process FASTP {
    tag "${sample_id}/${batch_id}"              // Show both IDs in log — much
    label 'process_medium'                      // more useful when debugging

    publishDir "${params.outdir}/fastp/${sample_id}", mode: 'copy',
        saveAs: { fn -> fn.endsWith('.fastq.gz') ? null : fn }
        // Publish only HTML and JSON (QC reports), not the trimmed reads.
        // Trimmed reads go straight to BWA_MEM2 via the channel — writing
        // them to disk and re-reading doubles I/O for no benefit.

    input:
    // Must match the 3-element tuple emitted by main.nf:
    // [sample_id, batch_id, [r1, r2]]
    tuple val(sample_id), val(batch_id), path(reads)

    output:
    tuple val(sample_id), val(batch_id), path("${sample_id}_${batch_id}_R1.trimmed.fastq.gz"),
                                          path("${sample_id}_${batch_id}_R2.trimmed.fastq.gz"), emit: reads
    path "${sample_id}_${batch_id}.fastp.html", emit: html
    path "${sample_id}_${batch_id}.fastp.json", emit: json

    script:
    def prefix = "${sample_id}_${batch_id}"    // Now actually useful — encodes both IDs
    """
    fastp \\
        --in1  ${reads[0]} \\
        --in2  ${reads[1]} \\
        --out1 ${prefix}_R1.trimmed.fastq.gz \\
        --out2 ${prefix}_R2.trimmed.fastq.gz \\
        --detect_adapter_for_pe \\
        --correction \\
        --low_complexity_filter \\
        --thread ${task.cpus} \\
        --qualified_quality_phred 20 \\
        --length_required 30 \\
        --html ${prefix}.fastp.html \\
        --json ${prefix}.fastp.json
    """

    stub:
    def prefix = "${sample_id}_${batch_id}"
    """
    touch ${prefix}_R1.trimmed.fastq.gz
    touch ${prefix}_R2.trimmed.fastq.gz
    touch ${prefix}.fastp.html
    touch ${prefix}.fastp.json
    """
}