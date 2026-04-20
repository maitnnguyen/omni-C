process BWA_MEM {
    tag "$sample_id"
    label 'process_high' // Alignment is CPU/RAM intensive

    publishDir "${params.outdir}/alignment", mode: 'copy',
        pattern: "*.{bam,bai}"

    input:
    // Matches the list output [R1, R2] from the improved FASTP module
    tuple val(sample_id), path(reads)
    path  index // Pass the genome index files as a path

    output:
    tuple val(sample_id), path("${sample_id}.bam"), emit: bam
    path "${sample_id}.bam.bai"                    , emit: bai

    script:
    // Omni-C recommendation: -5SP
    // -5: keep all parts of chimeric reads as separate alignments
    // -S: skip secondary alignment
    // -P: perform paired-end alignment (even if one end is chimeric)
    def prefix = "${sample_id}"
    """
    bwa mem \
        -5SP -M \
        -t ${task.cpus} \
        -R "@RG\\tID:${prefix}\\tSM:${prefix}\\tPL:ILLUMINA" \
        ${index}/genome.fa \
        ${reads[0]} \
        ${reads[1]} | \
        samtools sort -@ ${task.cpus} -o ${prefix}.bam -

    samtools index ${prefix}.bam
    """
}