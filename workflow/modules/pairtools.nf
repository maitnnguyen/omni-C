// modules/pairtools.nf

process PAIRTOOLS_PARSE {
    tag "${sample_id}/${batch_id}"
    label 'process_high'

    input:
    tuple val(sample_id), val(batch_id), path(bam)
    path  chrom_sizes
    // Genome-specific params passed as values from main.nf
    val   min_mapq
    val   walks_policy

    output:
    tuple val(sample_id), val(batch_id), path("${sample_id}_${batch_id}.pairs.gz"), emit: pairs
    path "${sample_id}_${batch_id}.parse.stats", emit: stats

    script:
    def prefix = "${sample_id}_${batch_id}"
    """
    pairtools parse \\
        --min-mapq    ${min_mapq} \\
        --walks-policy ${walks_policy} \\
        --max-inter-align-gap 30 \\
        --chroms-path ${chrom_sizes} \\
        --output-stats ${prefix}.parse.stats \\
        -o ${prefix}.pairs.gz \\
        ${bam}
    """

    stub:
    def prefix = "${sample_id}_${batch_id}"
    """
    touch ${prefix}.pairs.gz
    touch ${prefix}.parse.stats
    """
}


process PAIRTOOLS_MERGE_DEDUP {
    tag "${sample_id}"
    label 'process_high'
    // Sorting + dedup across merged batches is the most RAM/disk-intensive step
    // Label this 'process_high' and allocate accordingly in your config

    publishDir "${params.outdir}/pairs/${sample_id}", mode: 'copy'

    input:
    // groupTuple produces: [sample_id, [batch_ids...], [pairs_gz...]]
    tuple val(sample_id), val(batch_ids), path(batch_pairs)
    path  chrom_sizes

    output:
    tuple val(sample_id), path("${sample_id}.pairs.gz"),      emit: pairs
    path  "${sample_id}.dedup.stats",                          emit: stats

    script:
    """
    # Merge all per-batch pairs files for this sample into one stream,
    # sort across the merged set, then dedup.
    # Dedup MUST happen after merge — duplicate reads split across batches
    # from the same library will not be caught if you dedup per-batch.
    pairtools merge \\
        --nproc ${task.cpus} \\
        ${batch_pairs} \\
    | pairtools sort \\
        --nproc    ${task.cpus} \\
        --tmpdir   \${TMPDIR:-./} \\
    | pairtools dedup \\
        --output-stats ${sample_id}.dedup.stats \\
        -o ${sample_id}.pairs.gz
    """

    stub:
    """
    touch ${sample_id}.pairs.gz
    touch ${sample_id}.dedup.stats
    """
}