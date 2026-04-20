#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

// ── Imports — names must match process names defined in each module file ──────
include { FASTP                  } from './modules/fastp'
include { BWA_MEM2               } from './modules/bwa_mem2'
include { PAIRTOOLS_PARSE        } from './modules/pairtools'
include { PAIRTOOLS_MERGE_DEDUP  } from './modules/pairtools'
include { COOLER                 } from './modules/cooler'
include { JUICER                 } from './modules/juicer'
include { COOLTOOLS_SCALING      } from './modules/cooltools'

workflow {

    // ── Genome parameter resolution ───────────────────────────────────────────
    def genome_params = params.genomes[params.genome]
        ?: error("Unknown genome '${params.genome}'. " +
                 "Valid options: ${params.genomes.keySet().join(', ')}")

    // Index is a string prefix (not a path object) — bwa-mem2 resolves the
    // actual index files by appending .amb/.ann/.bwt.2bit.64 etc. itself.
    ch_index        = Channel.value(genome_params.bwa_index)

    // chrom_sizes is a real file — wrap in file() so Nextflow stages it.
    ch_chrom_sizes  = Channel.value(file(genome_params.chrom_sizes))

    // Genome-specific pairtools params — value channels so they broadcast
    // correctly to every parallel PAIRTOOLS_PARSE instance.
    ch_walks_policy = Channel.value(genome_params.walks_policy)
    ch_min_mapq     = Channel.value(genome_params.min_mapq)

    // ── 1. Input: parse metadata CSV ─────────────────────────────────────────
    // CSV columns: sample_id, batch_id, read1, read2
    // Tuple structure carried through the whole pipeline:
    //   per-batch processes : [sample_id, batch_id, ...]
    //   post-merge processes: [sample_id, ...]
    ch_reads = Channel
        .fromPath(params.input, checkIfExists: true)
        .splitCsv(header: true)
        .map { row ->
            tuple(
                row.sample_id,
                row.batch_id,
                file(row.read1, checkIfExists: true),
                file(row.read2, checkIfExists: true)
            )
        }

    // ── 2. Per-batch: trim ────────────────────────────────────────────────────
    // FASTP input : [sample_id, batch_id, r1, r2]
    // FASTP output: reads → [sample_id, batch_id, r1_trimmed, r2_trimmed]
    trimmed = FASTP(ch_reads)

    // ── 3. Per-batch: align ───────────────────────────────────────────────────
    // BWA_MEM2 input : [sample_id, batch_id, r1, r2], val index_prefix
    // BWA_MEM2 output: bam → [sample_id, batch_id, bam]
    bam = BWA_MEM2(trimmed.reads, ch_index)

    // ── 4. Per-batch: parse alignments to pairs ───────────────────────────────
    // PAIRTOOLS_PARSE input : [sample_id, batch_id, bam], chrom_sizes,
    //                         val min_mapq, val walks_policy
    // PAIRTOOLS_PARSE output: pairs → [sample_id, batch_id, pairs_gz]
    batch_pairs = PAIRTOOLS_PARSE(
        bam,
        ch_chrom_sizes,
        ch_min_mapq,
        ch_walks_policy
    )

    // ── 5. Per-sample: merge batches then dedup ───────────────────────────────
    // groupTuple(by: 0) collects all batches sharing the same sample_id:
    //   [sample_id, batch_id, pairs_gz] × N
    //   → [sample_id, [batch_id, ...], [pairs_gz, ...]]
    // Dedup MUST happen after merge — duplicates spanning batches from the
    // same library preparation will not be caught if you dedup per-batch.
    ch_to_merge = batch_pairs.pairs
        .groupTuple(by: 0)

    // PAIRTOOLS_MERGE_DEDUP output: pairs → [sample_id, pairs_gz]
    final_pairs = PAIRTOOLS_MERGE_DEDUP(ch_to_merge, ch_chrom_sizes)

    // ── 6. Matrix generation ──────────────────────────────────────────────────
    // Both tools consume the same final merged+deduped pairs per sample.
    // Access the named emit (.pairs) — without it you get the full output set.
    ch_to_cooler = final_pairs.pairs
        .groupTuple(by: 0)

    mcools = COOLER(final_pairs.pairs, ch_chrom_sizes)          // chrom_sizes 
    JUICER(final_pairs.pairs, ch_chrom_sizes)   // cooler reads genome from params

    // ── 7. QC ─────────────────────────────────────────────────────────────────
    // COOLTOOLS_SCALING takes the mcool output, not the pairs file.
    // chrom_sizes not needed by cooltools — it reads bin structure from mcool.
    COOLTOOLS_SCALING(mcools.mcool, ch_chrom_sizes)
}