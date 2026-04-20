process BWA_MEM2 {
    tag "${sample_id}/${batch_id}"
    label 'process_high'

    // scratch + index: only safe if your cluster config stages inputs to scratch.
    // If unsure, remove scratch and let the work dir handle it.
    scratch true

    input:
    // Match FASTP's 4-element output exactly
    tuple val(sample_id), val(batch_id), path(r1), path(r2)
    // Pass the resolved index prefix as a value string from params,
    // not as a path dir — avoids the runtime Groovy .find{} null problem.
    // In main.nf: Channel.value(genome_params.bwa_index)
    val index_prefix

    output:
    tuple val(sample_id), val(batch_id), path("${sample_id}_${batch_id}.bam"), emit: bam

    script:
    def prefix = "${sample_id}_${batch_id}"
    def rg     = "@RG\\tID:${prefix}\\tSM:${sample_id}\\tLB:${sample_id}\\tPU:${batch_id}\\tPL:ILLUMINA"
    """
    bwa-mem2 mem \\
        -5SP \\
        -t ${task.cpus} \\
        -R "${rg}" \\
        ${index_prefix} \\
        ${r1} \\
        ${r2} \\
    | samtools view \\
        -bS \\
        -@ ${task.cpus} \\
        -o ${prefix}.bam \\
        -
    """

    stub:
    def prefix = "${sample_id}_${batch_id}"
    """
    touch ${prefix}.bam
    """
}