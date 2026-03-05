#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// ============================================================================
//  WGS Metagenomics Pipeline — Oxford Nanopore Long-Read
//  ASD Biomarker-Focused Study
// ============================================================================

// ── Import modules ─────────────────────────────────────────────────────────
include { DORADO_BASECALL }        from './modules/basecalling'
include { NANOPLOT_QC }            from './modules/qc_nanoplot'
include { CHOPPER_FILTER }         from './modules/quality_filter'
include { MINIMAP2_HOST_REMOVAL }  from './modules/host_removal'
include { KRAKEN2_CLASSIFY }       from './modules/kraken2_classify'
include { BRACKEN_ABUNDANCE }      from './modules/bracken_abundance'
include { KRONA_PLOT }             from './modules/krona_visualization'
include { FASTQ_TO_FASTA }         from './modules/prodigal_genepred'
include { PRODIGAL_PREDICT }       from './modules/prodigal_genepred'
include { EGGNOG_ANNOTATE }        from './modules/eggnog_mapper'
include { EGGNOG_TO_PATHWAY_MATRIX } from './modules/eggnog_pathway_convert'
include { DIVERSITY_ANALYSIS }     from './modules/diversity_analysis'
include { LEFSE_ANALYSIS as LEFSE_TAXA }     from './modules/lefse_enrichment'
include { LEFSE_ANALYSIS as LEFSE_PATHWAYS } from './modules/lefse_enrichment'
include { MAASLIN2_ANALYSIS }      from './modules/maaslin2_correlation'
include { VISUALIZATION }          from './modules/visualization'

// ── Parameter validation ───────────────────────────────────────────────────
def validate_params() {
    if (!params.input_dir && !params.input_fastq) {
        error "ERROR: Provide either --input_dir (POD5/FAST5) or --input_fastq (FASTQ)"
    }
    if (!params.skip_basecalling && !params.input_dir) {
        error "ERROR: --input_dir is required when basecalling is enabled"
    }
    if (params.skip_basecalling && !params.input_fastq) {
        error "ERROR: --input_fastq is required when --skip_basecalling is true"
    }
    if (!params.kraken2_db) {
        error "ERROR: --kraken2_db is required (path to Kraken2 database directory)"
    }
    if (!params.skip_host_removal && !params.host_reference) {
        error "ERROR: --host_reference is required (GRCh38 FASTA/index) unless --skip_host_removal"
    }
    if (!params.skip_functional && !params.eggnog_db) {
        error "ERROR: --eggnog_db is required unless --skip_functional is true"
    }
}

// ── Log pipeline info ──────────────────────────────────────────────────────
log.info """
╔══════════════════════════════════════════════════════════════════════╗
║  WGS Metagenomics Pipeline — Nanopore Long-Read                    ║
║  ASD Biomarker Study                                               ║
╚══════════════════════════════════════════════════════════════════════╝

  Sample ID        : ${params.sample_id}
  Input dir        : ${params.input_dir ?: 'N/A'}
  Input FASTQ      : ${params.input_fastq ?: 'N/A'}
  Metadata         : ${params.metadata ?: 'N/A'}
  Output dir       : ${params.outdir}

  Skip basecalling : ${params.skip_basecalling}
  Skip host removal: ${params.skip_host_removal}
  Skip functional  : ${params.skip_functional}
  Skip diversity   : ${params.skip_diversity}
  Skip LEfSe       : ${params.skip_lefse}
  Skip MaAsLin2    : ${params.skip_maaslin2}

  Kraken2 DB       : ${params.kraken2_db ?: 'N/A'}
  Host reference   : ${params.host_reference ?: 'N/A'}
  eggNOG DB        : ${params.eggnog_db ?: 'N/A'}
  Max CPUs         : ${params.max_cpus}
  Max memory       : ${params.max_memory}
──────────────────────────────────────────────────────────────────────
""".stripIndent()

// ── Main workflow ──────────────────────────────────────────────────────────
workflow {

    validate_params()

    sample_id = params.sample_id

    // ── PHASE 1: Basecalling ───────────────────────────────────────────
    if (!params.skip_basecalling) {
        ch_pod5 = Channel.fromPath(params.input_dir, checkIfExists: true)
        DORADO_BASECALL(sample_id, ch_pod5)
        ch_raw_fastq = DORADO_BASECALL.out.fastq
    } else {
        ch_raw_fastq = Channel.fromPath(params.input_fastq, checkIfExists: true)
    }

    // ── PHASE 2: QC ───────────────────────────────────────────────────
    NANOPLOT_QC(sample_id, ch_raw_fastq)

    CHOPPER_FILTER(sample_id, ch_raw_fastq)
    ch_filtered = CHOPPER_FILTER.out.fastq

    // ── PHASE 3: Host removal ─────────────────────────────────────────
    if (!params.skip_host_removal) {
        ch_host_ref = Channel.fromPath(params.host_reference, checkIfExists: true)
        MINIMAP2_HOST_REMOVAL(sample_id, ch_filtered, ch_host_ref)
        ch_clean = MINIMAP2_HOST_REMOVAL.out.fastq
    } else {
        ch_clean = ch_filtered
    }

    // ── PHASE 4: Taxonomic classification ─────────────────────────────
    ch_kraken2_db = Channel.fromPath(params.kraken2_db, checkIfExists: true)

    KRAKEN2_CLASSIFY(sample_id, ch_clean, ch_kraken2_db)

    BRACKEN_ABUNDANCE(
        sample_id,
        KRAKEN2_CLASSIFY.out.report,
        ch_kraken2_db
    )

    KRONA_PLOT(sample_id, KRAKEN2_CLASSIFY.out.report)

    ch_abundance = BRACKEN_ABUNDANCE.out.abundance

    // ── PHASE 5: Functional profiling ─────────────────────────────────
    if (!params.skip_functional) {
        FASTQ_TO_FASTA(sample_id, ch_clean)
        PRODIGAL_PREDICT(sample_id, FASTQ_TO_FASTA.out.fasta)

        ch_eggnog_db = Channel.fromPath(params.eggnog_db, checkIfExists: true)
        EGGNOG_ANNOTATE(
            sample_id,
            PRODIGAL_PREDICT.out.proteins,
            ch_eggnog_db
        )

        // Convert raw eggNOG annotations → KEGG pathway abundance table
        // (required for LEfSe pathways — eggNOG and Bracken have different formats)
        EGGNOG_TO_PATHWAY_MATRIX(sample_id, EGGNOG_ANNOTATE.out.annotations)
        ch_pathway_abundance = EGGNOG_TO_PATHWAY_MATRIX.out.pathway_abundance
    }

    // ── PHASE 6: Diversity analysis ───────────────────────────────────
    if (!params.skip_diversity && params.metadata) {
        ch_metadata = Channel.fromPath(params.metadata, checkIfExists: true)
        DIVERSITY_ANALYSIS(sample_id, ch_abundance, ch_metadata)
        ch_diversity = DIVERSITY_ANALYSIS.out.results
    }

    // ── PHASE 7: LEfSe differential enrichment ───────────────────────
    if (!params.skip_lefse && params.metadata) {
        ch_metadata_lefse = Channel.fromPath(params.metadata, checkIfExists: true)

        LEFSE_TAXA(sample_id, ch_abundance, ch_metadata_lefse, 'taxa')

        if (!params.skip_functional) {
            LEFSE_PATHWAYS(sample_id, ch_pathway_abundance, ch_metadata_lefse, 'pathways')
        }
    }

    // ── PHASE 8: MaAsLin2 biomarker-microbiome correlation ───────────
    if (!params.skip_maaslin2 && params.metadata) {
        ch_metadata_maaslin = Channel.fromPath(params.metadata, checkIfExists: true)
        MAASLIN2_ANALYSIS(sample_id, ch_abundance, ch_metadata_maaslin)
        ch_maaslin2 = MAASLIN2_ANALYSIS.out.results
    }

    // ── PHASE 9: Visualization ───────────────────────────────────────
    if (!params.skip_visualization && params.metadata) {
        ch_metadata_viz = Channel.fromPath(params.metadata, checkIfExists: true)

        // Collect results from upstream steps (use dummy if skipped)
        div_results  = (!params.skip_diversity  && params.metadata) ? ch_diversity  : Channel.fromPath("${projectDir}/assets/PLACEHOLDER")
        maa_results  = (!params.skip_maaslin2   && params.metadata) ? ch_maaslin2   : Channel.fromPath("${projectDir}/assets/PLACEHOLDER")

        VISUALIZATION(
            sample_id,
            ch_abundance,
            ch_metadata_viz,
            div_results,
            maa_results
        )
    }
}

// ── Completion handler ─────────────────────────────────────────────────────
workflow.onComplete {
    log.info """
    ──────────────────────────────────────────────────────────────────
    Pipeline completed at : ${workflow.complete}
    Duration              : ${workflow.duration}
    Success               : ${workflow.success}
    Output directory      : ${params.outdir}
    ──────────────────────────────────────────────────────────────────
    """.stripIndent()
}
