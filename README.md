# WGS Metagenomics Pipeline — Oxford Nanopore Long-Read

Whole-genome shotgun metagenomics pipeline for the ASD Biomarker Study.
Built with Nextflow DSL2 using modular per-process Docker containers.

## Directory Structure

```
nextflow-codebase/
├── main.nf                          # Main pipeline entry point
├── nextflow.config                  # Pipeline configuration & defaults
├── conf/
│   └── test.config                  # Test profile (minimal run)
├── modules/
│   ├── basecalling.nf               # Dorado basecaller
│   ├── qc_nanoplot.nf               # NanoPlot quality control
│   ├── quality_filter.nf            # Chopper read filtering
│   ├── host_removal.nf              # Minimap2 host DNA removal
│   ├── kraken2_classify.nf          # Kraken2 taxonomic classification
│   ├── bracken_abundance.nf         # Bracken species abundance
│   ├── krona_visualization.nf       # Krona interactive taxonomy plot
│   ├── prodigal_genepred.nf         # Prodigal gene prediction
│   ├── eggnog_mapper.nf             # eggNOG-mapper functional annotation
│   ├── diversity_analysis.nf        # Alpha/beta diversity + ANCOM-BC
│   ├── lefse_enrichment.nf          # LEfSe differential enrichment
│   ├── maaslin2_correlation.nf      # MaAsLin2 biomarker correlation
│   └── visualization.nf             # Publication-quality figures
├── bin/
│   ├── diversity_analysis.R         # R script for diversity analysis
│   ├── maaslin2_analysis.R          # R script for MaAsLin2
│   ├── visualization.R              # R script for figure generation
│   └── format_lefse_input.py        # Python script to format LEfSe input
├── docker/
│   └── r-metagenomics/Dockerfile    # Custom R image (only image you build)
└── assets/
    ├── sample_metadata_template.csv # Metadata template
    └── PLACEHOLDER                  # Placeholder for optional inputs
```

## Prerequisites

### 1. Install Nextflow

```bash
# Requires Java 11+
curl -s https://get.nextflow.io | bash
sudo mv nextflow /usr/local/bin/
```

### 2. Install Docker

```bash
# macOS
brew install --cask docker

# Linux (Ubuntu)
sudo apt-get update
sudo apt-get install docker.io
sudo usermod -aG docker $USER
```

### 3. Build the custom R Docker image

This is the **only image you need to build**. All other tools use public images.

```bash
cd docker/r-metagenomics
docker build -t bio-r-metagenomics:1.0.0 .
cd ../..
```

### 4. Download required databases

```bash
# ── Kraken2 + Bracken database (Standard + Protozoa & Fungi, ~16 GB) ──
mkdir -p databases/kraken2_pluspf
cd databases/kraken2_pluspf
wget https://genome-idx.s3.amazonaws.com/kraken/k2_pluspf_16gb_20231009.tar.gz
tar xzf k2_pluspf_16gb_20231009.tar.gz
cd ../..

# ── GRCh38 human reference genome (for host removal) ──
mkdir -p databases/grch38
cd databases/grch38
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/seqs_for_alignment_pipelines.ucsc_ids/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz
gunzip GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz
mv GCA_000001405.15_GRCh38_no_alt_analysis_set.fna GRCh38.fa
cd ../..

# ── eggNOG database (~45 GB, only if running functional profiling) ──
mkdir -p databases/eggnog
cd databases/eggnog
download_eggnog_data.py --data_dir .
cd ../..
```

## Running the Pipeline

### Full pipeline (from raw POD5 files)

```bash
nextflow run main.nf \
    --input_dir       /path/to/pod5_files/ \
    --metadata        /path/to/sample_metadata.csv \
    --host_reference  /path/to/databases/grch38/GRCh38.fa \
    --kraken2_db      /path/to/databases/kraken2_pluspf/ \
    --eggnog_db       /path/to/databases/eggnog/ \
    --sample_id       SAMPLE001 \
    --outdir          ./results
```

### From pre-basecalled FASTQ (skip Dorado)

```bash
nextflow run main.nf \
    --skip_basecalling \
    --input_fastq     /path/to/reads.fastq \
    --metadata        /path/to/sample_metadata.csv \
    --host_reference  /path/to/databases/grch38/GRCh38.fa \
    --kraken2_db      /path/to/databases/kraken2_pluspf/ \
    --eggnog_db       /path/to/databases/eggnog/ \
    --sample_id       SAMPLE001 \
    --outdir          ./results
```

### Taxonomy-only run (skip functional + downstream)

```bash
nextflow run main.nf \
    --skip_basecalling \
    --input_fastq     /path/to/reads.fastq \
    --host_reference  /path/to/databases/grch38/GRCh38.fa \
    --kraken2_db      /path/to/databases/kraken2_pluspf/ \
    --skip_functional \
    --skip_diversity \
    --skip_lefse \
    --skip_maaslin2 \
    --skip_visualization \
    --sample_id       SAMPLE001 \
    --outdir          ./results
```

### Dry run (see what would execute without running)

```bash
nextflow run main.nf -preview \
    --skip_basecalling \
    --input_fastq     /path/to/reads.fastq \
    --host_reference  /path/to/databases/grch38/GRCh38.fa \
    --kraken2_db      /path/to/databases/kraken2_pluspf/ \
    --skip_functional \
    --sample_id       TEST
```

### Resume a failed/interrupted run

```bash
nextflow run main.nf -resume [same params as before]
```

## Parameters

| Parameter | Default | Description |
|---|---|---|
| `--input_dir` | null | Directory with POD5/FAST5 raw signal files |
| `--input_fastq` | null | Pre-basecalled FASTQ file (use with `--skip_basecalling`) |
| `--metadata` | null | Sample metadata CSV (see template in `assets/`) |
| `--outdir` | `./results` | Output directory |
| `--sample_id` | `sample` | Sample identifier used in file naming |
| `--skip_basecalling` | false | Skip Dorado basecalling (provide `--input_fastq`) |
| `--dorado_model` | `sup` | Dorado model: `sup`, `hac`, or `fast` |
| `--min_quality` | 10 | Chopper minimum quality score |
| `--min_length` | 500 | Chopper minimum read length (bp) |
| `--host_reference` | null | GRCh38 FASTA for host removal |
| `--skip_host_removal` | false | Skip host decontamination |
| `--kraken2_db` | null | Kraken2 database directory |
| `--bracken_read_length` | 1500 | Bracken average read length |
| `--bracken_level` | `S` | Bracken level: S(pecies), G(enus), F(amily) |
| `--bracken_threshold` | 10 | Minimum reads for Bracken inclusion |
| `--skip_functional` | false | Skip Prodigal + eggNOG-mapper |
| `--eggnog_db` | null | eggNOG database directory |
| `--skip_diversity` | false | Skip diversity analysis |
| `--fdr_alpha` | 0.05 | FDR cutoff for diversity tests |
| `--fdr_ancombc` | 0.2 | FDR cutoff for ANCOM-BC |
| `--skip_lefse` | false | Skip LEfSe analysis |
| `--lda_threshold` | 2.0 | LEfSe LDA score threshold |
| `--skip_maaslin2` | false | Skip MaAsLin2 correlation |
| `--fixed_effects` | `diagnosis,age,sex` | MaAsLin2 fixed effects |
| `--reference_level` | `neurotypical` | MaAsLin2 reference group |
| `--max_cpus` | 16 | Maximum CPUs per process |
| `--max_memory` | `64.GB` | Maximum memory per process |

## Output Structure

```
results/
├── 01_qc/
│   └── nanoplot/              # NanoPlot QC reports and plots
├── 02_host_removal/           # Host removal statistics
├── 03_taxonomy/
│   ├── kraken2/               # Kraken2 classification output + report
│   ├── bracken/               # Bracken species-level abundance tables
│   └── krona/                 # Interactive Krona HTML taxonomy plot
├── 04_functional/
│   ├── prodigal/              # Predicted proteins and genes
│   └── eggnog/                # KEGG, COG, GO annotations
├── 05_diversity/              # Alpha/beta diversity, ANCOM-BC results
├── 06_lefse/                  # LEfSe differential enrichment results
├── 07_maaslin2/               # Biomarker-microbiome correlations
├── 08_visualization/          # Publication-quality figures
└── pipeline_info/
    ├── timeline.html          # Execution timeline
    ├── report.html            # Resource usage report
    └── dag.html               # Pipeline DAG visualization
```

## Docker Images Used

| Process | Image | Source | Build? |
|---|---|---|---|
| Dorado | `ontresearch/dorado:latest` | Docker Hub | No |
| NanoPlot | `staphb/nanoplot:latest` | Docker Hub | No |
| Chopper | `quay.io/biocontainers/chopper:0.9.0` | BioContainers | No |
| Minimap2 + samtools | `staphb/minimap2:latest` | Docker Hub | No |
| Kraken2 | `staphb/kraken2:latest` | Docker Hub | No |
| Bracken | `staphb/bracken:latest` | Docker Hub | No |
| Krona | `nanozoo/krona:2.7.1` | Docker Hub | No |
| Prodigal | `biocontainers/prodigal:v2.6.3` | BioContainers | No |
| eggNOG-mapper | `nanozoo/eggnog-mapper:2.1.12` | Docker Hub | No |
| LEfSe | `biobakery/lefse:latest` | Docker Hub | No |
| R metagenomics | `bio-r-metagenomics:1.0.0` | Local build | **Yes** |

## Metadata File Format

See `assets/sample_metadata_template.csv`. Required columns:

- `sample_id` — unique sample identifier
- `group` — `ASD` or `neurotypical`
- `age` — age in years (5-9)
- `sex` — `M` or `F`
- Biomarker columns: `vitamin_d`, `homocysteine`, `iron`, `tsh`, `total_cholesterol`, `triglycerides`, `hdl`, `ldl`, `calprotectin`, `lactoferrin`
- Clinical columns: `gi_symptoms`, `asd_severity`, `diet_type`

## Multi-Sample Cohort Analysis

This pipeline processes **one sample at a time**. For the full 204-sample cohort analysis:

1. Run the pipeline for each sample individually (basecalling through Bracken)
2. Merge all Bracken output tables into a single abundance matrix
3. Re-run the diversity, LEfSe, MaAsLin2, and visualization steps with the merged matrix and full metadata

To process multiple samples in a loop:

```bash
for FASTQ in /path/to/fastqs/*.fastq; do
    SAMPLE=$(basename "$FASTQ" .fastq)
    nextflow run main.nf \
        --skip_basecalling \
        --input_fastq "$FASTQ" \
        --host_reference /path/to/GRCh38.fa \
        --kraken2_db /path/to/kraken2_pluspf/ \
        --skip_functional --skip_diversity --skip_lefse --skip_maaslin2 --skip_visualization \
        --sample_id "$SAMPLE" \
        --outdir "./results/${SAMPLE}" \
        -resume
done
```
