# Omni-C sequencing pipeline

**Omni-C · Contact Matrix · SNV/SV Calling · Haplotype Phasing**

A Nextflow DSL2 pipeline for processing Omni-C sequencing data from raw FASTQ through contact matrix generation, variant calling (SNV/indel + SV), and chromosome-scale haplotype phasing.

**Author:** Mai TN Nguyen · [@maitnnguyen](https://github.com/maitnnguyen) · ntnmai303@gmail.com  
**Version:** 3.0.0 · **Nextflow:** ≥ 23.04.0 · **Sample type:** Human PBMC / diploid

---

## Overview

Omni-C uses a sequence-independent endonuclease (DNase-based) to fragment chromatin, producing uniform genome-wide coverage without restriction site bias. This makes a single Omni-C library suitable for multiple analyses simultaneously.

This pipeline takes raw paired-end FASTQ files and produces:

| Output | Tool | Location |
|---|---|---|
| ICE-balanced contact matrix (`.mcool`) | cooler | `results/contacts/` |
| Filtered SNV + indel VCF | GATK4 | `results/variants/filtered/` |
| Merged SV VCF | LUMPY + DELLY + SURVIVOR | `results/sv/merged/` |
| Phased haplotype blocks + VCF | HapCUT2 `--hic` | `results/phasing/hapcut2/` |
| Final combined VCF (SNV + SV + phase tags) | bcftools | `results/final/` |
| Aggregate QC report | MultiQC | `results/qc/multiqc/` |

---

## Pipeline steps

```
Raw FASTQ (R1 + R2)
        │
        ├── [1] fastp                   Read QC + auto detect adapter trimming
        │
        ├── [2] BWA-MEM2 (-5SP)         Align in Omni-C chimeric read mode
        │
        ├── [3] samtools                Sort + index BAM
        │
        ├── [4] pairtools               Parse ligation junctions → dedup → stats
        │
        ├── [5] cooler                  Bin contacts → .mcool (ICE balanced)
        │
        ├── [6] GATK MarkDuplicates     Mark PCR/optical duplicates
        │
        ├── [7] GATK BQSR (optional)    Base quality recalibration
        │
        ├── [8–10] GATK HaplotypeCaller SNV/indel calling (GVCF → filter)
        │
        ├── [11–12] LUMPY + DELLY       SV calling (split reads + discordant pairs)
        │
        ├── [13] SURVIVOR               Merge SV calls across callers
        │
        ├── [14] HapCUT2 --hic          Extract fragment matrix + phase variants
        │
        ├── [15] bcftools               Merge phased SNV + SV into final VCF
        │
        └── [16] MultiQC               Aggregate QC
```

---

## Reference genome options

The pipeline works with any indexed reference genome. Two builds are commonly used for human data:

### GRCh38 (hg38)

The standard GATK bundle reference, widely used and well-supported by known-sites VCFs for BQSR.

```bash
# FASTA
wget https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/Homo_sapiens_assembly38.fasta

# BWA-MEM2 index
bwa-mem2 index Homo_sapiens_assembly38.fasta

# Known sites for BQSR
wget https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/Homo_sapiens_assembly38.dbsnp138.vcf
wget https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz
```

### T2T-CHM13v2.0 (hs1) — recommended for Omni-C

The complete telomere-to-telomere human reference (T2T Consortium, 2022) resolves ~200 Mb of sequence absent from GRCh38, including centromeres, pericentromeric regions, and the short arms of acrocentric chromosomes. For Omni-C, T2T is particularly valuable because Omni-C's uniform (non-restriction-biased) coverage reaches these previously inaccessible regions. Benefits include:

- Eliminates reference gaps that cause spurious contacts and mapping artefacts in heterochromatic regions
- Improves SNV calling accuracy in segmental duplications and repeat-rich regions
- Reduces mapping bias at centromeres and telomeres, which are captured by Omni-C but not by restriction-enzyme Hi-C
- Pre-built bwa-mem2 index available from UCSC — no need to build from scratch

```bash
# Download FASTA (chrN naming, recommended)
wget https://hgdownload.soe.ucsc.edu/hubs/GCA/009/914/755/GCA_009914755.4/GCA_009914755.4.chrNames.fa.gz
gunzip GCA_009914755.4.chrNames.fa.gz

# Download pre-built bwa-mem2 index (saves ~1 h of compute)
wget https://hgdownload.soe.ucsc.edu/hubs/GCA/009/914/755/GCA_009914755.4/bwa-mem2/GCA_009914755.4.chrNames.fa.gz.bwt.2bit.64
wget https://hgdownload.soe.ucsc.edu/hubs/GCA/009/914/755/GCA_009914755.4/bwa-mem2/GCA_009914755.4.chrNames.fa.gz.pac
wget https://hgdownload.soe.ucsc.edu/hubs/GCA/009/914/755/GCA_009914755.4/bwa-mem2/GCA_009914755.4.chrNames.fa.gz.0123
wget https://hgdownload.soe.ucsc.edu/hubs/GCA/009/914/755/GCA_009914755.4/bwa-mem2/GCA_009914755.4.chrNames.fa.gz.amb
wget https://hgdownload.soe.ucsc.edu/hubs/GCA/009/914/755/GCA_009914755.4/bwa-mem2/GCA_009914755.4.chrNames.fa.gz.ann
```

> **Note on BQSR with T2T:** Known-sites VCFs for T2T are available from the 1000 Genomes Project T2T release on AnVIL (`AnVIL_T2T_CHRY` workspace). If these are unavailable, BQSR can be skipped — the pipeline handles this gracefully when `--known_snps` and `--known_indels` are not provided.

---

## Installation

### 1. Nextflow

```bash
curl -s https://get.nextflow.io | bash
mv nextflow ~/bin/
nextflow -version   # requires >= 23.04.0
```

### 2. Singularity / Apptainer

Confirm available on your HPC:
```bash
singularity --version   # or: apptainer --version
```

### 3. Clone this repository

```bash
git clone https://github.com/maitnnguyen/omni-C.git
cd omni-C
```

### 4. Pull Singularity SIF images

All tools are containerized. Pull the SIF files manually before running. GATK and samtools are loaded via HPC modules (see [Container setup](#container-setup)).

```bash
cd containers/sif/

# Alignment + contacts
singularity pull --disable-cache bwa_mem2_samtools.sif    docker://quay.io/biocontainers/mulled-v2-e5d375990341c5aef3c9aff74f96f66f65375ef6:c5b8c4b7735290369693e2b63cfc1ea0732fde07-0
singularity pull --disable-cache pairtools_1.1.0.sif      docker://bskubi/pairtools:1.1.0
singularity pull --disable-cache cooler_0.9.3.sif         docker://quay.io/biocontainers/cooler:0.9.3--pyhdfd78af_0

# QC
singularity pull --disable-cache fastqc_0.12.1.sif        docker://quay.io/biocontainers/fastqc:0.12.1--hdfd78af_0
singularity pull --disable-cache fastp_0.23.4.sif         docker://quay.io/biocontainers/fastp:0.23.4--hadf994f_2
singularity pull --disable-cache multiqc_1.19.sif         docker://quay.io/biocontainers/multiqc:1.19--pyhdfd78af_0

# SV calling
singularity pull --disable-cache lumpy_0.3.1.sif          docker://quay.io/biocontainers/lumpy-sv:0.3.1-3
singularity pull --disable-cache delly_1.1.6.sif          docker://quay.io/biocontainers/delly:1.1.6--ha41ced6_1
singularity pull --disable-cache survivor_1.0.7.sif       docker://quay.io/biocontainers/survivor:1.0.7--hec16e2b_3

# Phasing + merging
singularity pull --disable-cache hapcut2_1.3.4.sif        docker://quay.io/biocontainers/hapcut2:1.3.4-2
singularity pull --disable-cache bcftools_1.18.sif        docker://quay.io/biocontainers/bcftools:1.18--h8b25389_0
```

---

## Usage

### Samplesheet

Create a CSV file with one row per sample:

```csv
sample, [batch],fastq_1,fastq_2
SAMPLE1,[batch],/data/SAMPLE1_R1.fastq.gz,/data/SAMPLE1_R2.fastq.gz
SAMPLE2,[batch],/data/SAMPLE2_R1.fastq.gz,/data/SAMPLE2_R2.fastq.gz
```

### Run — GRCh38

```bash
nextflow run main.nf \
    -profile slurm,singularity \
    --input    samplesheet.csv \
    --genome    hg38 \ 
    --known_snps   /ref/GRCh38/Homo_sapiens_assembly38.dbsnp138.vcf \
    --known_indels /ref/GRCh38/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz \
    --genome   GRCh38 \
    --outdir   results/
```

### Run — T2T-CHM13v2.0

```bash
nextflow run main.nf \
    -profile slurm,singularity \
    --input      samplesheet.csv \
    --genome     t2t \
    --bwa_index  /ref/T2T/bwa-mem2/GCA_009914755.4.chrNames.fa \
    --genome     T2T-CHM13v2 \
    --outdir     results_t2t/
```

> When running with T2T and no known-sites VCFs, the pipeline skips BQSR automatically and logs a warning. This is expected behaviour.

### Run — local Docker (testing)

```bash
nextflow run main.nf \
    -profile docker \
    --input  samplesheet.csv \
    --genome  hg38 \
    --outdir results/
```

### Resume after interruption

```bash
nextflow run main.nf -profile slurm,singularity --input ... -resume
```

---

## Key parameters

| Parameter | Default | Description |
|---|---|---|
| `--input` | required | Samplesheet CSV |
| `--genome`| required | Reference genome build |
| `--bwa_index` | auto-built | Path to BWA-MEM2 index prefix (skip build if pre-downloaded) |
| `--known_snps` | null | dbSNP VCF for BQSR (optional) |
| `--known_indels` | null | Mills+1000G indels VCF for BQSR (optional) |
| `--sif_dir` | `./containers/sif` | Directory containing pre-pulled SIF files |
| `--module_gatk` | `GATK/4.4.0.0` | HPC module name for GATK |
| `--module_samtools` | `samtools/1.18` | HPC module name for samtools |
| `--bwa_extra` | `-5SP` | BWA-MEM2 flags — **do not change** for Omni-C |
| `--mapq_filter` | `30` | Minimum MAPQ for valid contact pairs, 40 for t2t2 |
| `--cooler_resolutions` | `1000,...,1000000` | Bin sizes for `.mcool` output |
| `--sv_merge_dist` | `1000` | Max bp distance to merge SV calls across callers |
| `--outdir` | `results` | Output directory |
| `--max_cpus` | `16` | Max CPUs per process |
| `--max_memory` | `128.GB` | Max memory per process |

---

## Container setup

### Tools loaded via HPC module (no SIF required)

| Tool | Module | Processes |
|---|---|---|
| GATK 4.4.0.0 | `GATK/4.4.0.0` | MarkDuplicates, BQSR, HaplotypeCaller, GenotypeGVCFs, VariantFiltration |
| samtools 1.18 | `samtools/1.18` | sort, index, stats, faidx, chromsizes |

Adjust module names in `nextflow.config` to match your HPC:
```groovy
params.module_gatk     = "GATK/4.4.0.0"
params.module_samtools = "samtools/1.18"
```

### Tools using SIF containers

| Tool | Version | SIF filename | Image |
|---|---|---|---|
| BWA-MEM2 + samtools | 2.2.1 + 1.13 | `bwa_mem2_samtools.sif` | mulled biocontainer |
| pairtools | 1.1.0 | `pairtools_1.1.0.sif` | `bskubi/pairtools:1.1.0` |
| cooler | 0.9.3 | `cooler_0.9.3.sif` | quay.io/biocontainers |
| FastQC | 0.12.1 | `fastqc_0.12.1.sif` | quay.io/biocontainers |
| fastp | 0.23.4 | `fastp_0.23.4.sif` | quay.io/biocontainers |
| MultiQC | 1.19 | `multiqc_1.19.sif` | quay.io/biocontainers |
| LUMPY | 0.3.1 | `lumpy_0.3.1.sif` | quay.io/biocontainers |
| DELLY | 1.1.6 | `delly_1.1.6.sif` | quay.io/biocontainers |
| SURVIVOR | 1.0.7 | `survivor_1.0.7.sif` | quay.io/biocontainers |
| HapCUT2 | 1.3.4 | `hapcut2_1.3.4.sif` | quay.io/biocontainers |
| bcftools | 1.18 | `bcftools_1.18.sif` | quay.io/biocontainers |

> **Why a mulled container for BWA-MEM2?** The alignment step pipes `bwa-mem2` output directly into `samtools view`. Both tools must be available in the same container environment. The mulled biocontainer `mulled-v2-e5d375990341c5aef3c9aff74f96f66f65375ef6` bundles both tools together. HPC modules are not accessible inside Singularity containers, which is why a combined container is needed here specifically.

---

## Output structure

```
results/
├── qc/
│   ├── fastqc/              Per-sample FastQC reports
│   ├── fastp/               Trimming stats
│   ├── samtools_stats/      Alignment statistics
│   └── multiqc/             ⭐ Open multiqc_report.html for full QC overview
├── contacts/
│   └── SAMPLE/
│       ├── *.dedup.stats.txt     pairtools QC: cis/trans ratio, pair types
│       ├── *.library.stats.txt   Full library statistics
│       └── *.mcool               ⭐ Multi-resolution contact matrix (open in HiGlass/Juicebox)
├── variants/
│   ├── markdup/             MarkDuplicates metrics
│   ├── gvcf/                Per-sample GVCFs
│   ├── genotyped/           Joint-genotyped VCFs
│   └── filtered/            ⭐ Hard-filtered SNV + indel VCFs
├── sv/
│   ├── lumpy/               LUMPY raw calls
│   ├── delly/               DELLY raw calls
│   └── merged/              ⭐ SURVIVOR-merged SV VCF + per-type stats
├── phasing/
│   └── hapcut2/
│       ├── *.haplotype_blocks    HapCUT2 phase blocks (text format)
│       ├── *.phased.vcf.gz       Phased VCF with PS/HP tags
│       └── *.phase_summary.txt   Block count, variants phased, largest blocks
├── final/
│   └── SAMPLE/
│       ├── *.final.vcf.gz        ⭐ Phased SNVs + SVs in one VCF
│       └── *.summary.txt         Variant counts by type
└── pipeline_info/
    ├── timeline.html        Execution timeline
    ├── report.html          Resource usage report
    ├── trace.txt            Per-task trace
    └── dag.html             Pipeline DAG
```

---

## Notes on Omni-C alignment flags

BWA-MEM2 is run with `-5SP` (equivalent to the `-5 -S -P` combination). These flags are **required** for Omni-C and must not be changed:

| Flag | Effect | Why it matters for Omni-C |
|---|---|---|
| `-5` | Report the 5-prime segment as primary for chimeric reads | Ensures pairtools correctly identifies ligation junctions |
| `-S` | Skip mate rescue | The two read ends are NOT true mates — they come from different genomic loci |
| `-P` | Skip pairing | Treats each end independently, required for proximity ligation data |

---

## Haplotype phasing with HapCUT2

This pipeline uses HapCUT2 with the `--hic` flag, which is specifically designed to exploit the long-range linkage information in Omni-C (and Hi-C) data.

Each chimeric Omni-C read pair connects two genomic loci that can be megabases apart on the same chromosome. `extractHAIRS --hic` converts these pairs into a fragment matrix that encodes haplotype co-occurrence across the full distance, rather than just within a single short read. The result is phase blocks orders of magnitude larger than standard short-read phasing.

**Expected output quality (human, ~200M valid pairs):**

| Metric | Typical value |
|---|---|
| Phase block N50 | 20–60 Mb |
| % heterozygous SNVs phased | 85–95% |
| Switch error rate | < 0.5% |

---

## Citation

If you use this pipeline, please cite the underlying tools:

- **pairtools**: Open2C et al., *PLOS Computational Biology* (2023)
- **cooler**: Abdennur & Mirny, *Bioinformatics* (2020)
- **GATK4**: McKenna et al., *Genome Research* (2010)
- **LUMPY**: Layer et al., *Genome Biology* (2014)
- **DELLY**: Rausch et al., *Bioinformatics* (2012)
- **SURVIVOR**: Jeffares et al., *Nature Communications* (2017)
- **HapCUT2**: Edge et al., *Genome Research* (2017)
- **T2T-CHM13**: Nurk et al., *Science* (2022)

---

## License

MIT