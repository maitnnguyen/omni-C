// modules/juicer.nf

process JUICER {
    tag "${sample_id}"
    label 'process_high'

    publishDir "${params.outdir}/hic/${sample_id}", mode: 'copy'

    input:
    tuple val(sample_id), path(pairs)
    path  chrom_sizes

    output:
    tuple val(sample_id), path("${sample_id}.hic"), emit: hic

    script:
    def prefix   = "${sample_id}"
    def hic_res  = params.hic_resolutions ?: "1000,5000,10000,25000,50000,100000,250000,500000,1000000"
    // Give JVM 80% of allocated memory; keep 20% headroom for JVM overhead.
    // task.memory is a nextflow Memory object; toMega() gives integer MB.
    def jvm_mem  = task.memory ? "-Xmx${(task.memory.toMega() * 0.8).intValue()}m" : "-Xmx26000m"
    // juicer_tools jar path comes from params — never hardcode tool paths.
    // Set params.juicer_jar in nextflow.config or pass --juicer_jar on CLI.
    def jar      = params.juicer_jar
    """
    # pairtools .pairs.gz is not directly readable by juicer_tools pre.
    # Reformat to the short-format pairs juicer expects:
    # readID str1 pos1 str2 pos2 frag1 frag2
    # We use dummy fragment IDs (0/1) since we are not restriction-site based.
    pairtools select 'pair_type=="UU" or pair_type=="RU" or pair_type=="UR"' \\
        ${pairs} \\
    | pairtools flip \\
    | pairtools header dump --header - \\
    | paste - - \\
    | awk 'BEGIN{OFS="\\t"} !/^#/ {
        print \$1, \$2, \$3, \$5, \$4, \$6, 0, 1
      }' \\
    | sort -k2,2d -k6,6d \\
    | gzip > ${prefix}.juicer.txt.gz

    java ${jvm_mem} -jar ${jar} pre \\
        -r ${hic_res} \\
        ${prefix}.juicer.txt.gz \\
        ${prefix}.hic \\
        ${chrom_sizes}

    rm ${prefix}.juicer.txt.gz
    """

    stub:
    """
    touch ${sample_id}.hic
    """
}