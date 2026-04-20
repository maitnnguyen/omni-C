process PAIRTOOLS_MERGE_DEDUP {
    tag "$sample_id"
    label 'process_high'
    
    // Crucial for heavy data: use local SSD for sorting/merging
    scratch true 

    publishDir "${params.outdir}/merged_pairs", mode: 'copy'

    input:
    tuple val(sample_id), path(batch_pairs)
    path  chrom_sizes

    output:
    tuple val(sample_id), path("${sample_id}.final.pairs.gz"), emit: pairs
    path "${sample_id}.merged.pairstat.txt"                 , emit: stats

    script:
    def prefix = "${sample_id}"
    """
    # 1. Merge all batch files into one stream
    # 2. Sort the merged stream (required before dedup)
    # 3. Dedup the sorted stream
    
    pairtools merge \
        --nproc ${task.cpus} \
        --chroms-path ${chrom_sizes} \
        ${batch_pairs} | \
    pairtools sort \
        --nproc ${task.cpus} \
        --tmpdir ./ \
        --memory ${task.memory.toGiga()}G | \
    pairtools dedup \
        --output-stats ${prefix}.merged.pairstat.txt \
        --output ${prefix}.final.pairs.gz
    """
}