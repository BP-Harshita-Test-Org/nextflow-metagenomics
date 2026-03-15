#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// ============================================================================
//  WGS Metagenomics Pipeline — Oxford Nanopore Long-Read
//  ASD Biomarker-Focused Study (Multi-Sample)
// ============================================================================

// ── Import modules ─────────────────────────────────────────────────────────
include { NANOPLOT_QC }            from './modules/qc_nanoplot'
include { CHOPPER_FILTER }         from './modules/quality_filter'
include { MINIMAP2_HOST_REMOVAL }  from './modules/host_removal'
include { KRAKEN2_CLASSIFY }       from './modules/kraken2_classify'
include { BRACKEN_ABUNDANCE }      from './modules/bracken_abundance'
include { KRONA_PLOT }             from './modules/krona_visualization'
include { METAFLYE_ASSEMBLE }      from './modules/metaflye_assembly'
include { QUAST_ASSESSMENT }       from './modules/metaflye_assembly'
include { SEQKIT_FILTER }          from './modules/metaflye_assembly'
include { FASTQ_TO_FASTA }         from './modules/prodigal_genepred'
include { PRODIGAL_PREDICT }       from './modules/prodigal_genepred'
include { EGGNOG_ANNOTATE }        from './modules/eggnog_mapper'
include { EGGNOG_TO_PATHWAY_MATRIX } from './modules/eggnog_pathway_convert'
include { MERGE_ABUNDANCE as MERGE_TAXA }       from './modules/merge_abundance'
include { MERGE_ABUNDANCE as MERGE_PATHWAYS }   from './modules/merge_abundance'
include { DIVERSITY_ANALYSIS }     from './modules/diversity_analysis'
include { LEFSE_ANALYSIS as LEFSE_TAXA }     from './modules/lefse_enrichment'
include { LEFSE_ANALYSIS as LEFSE_PATHWAYS } from './modules/lefse_enrichment'
include { MAASLIN2_ANALYSIS }      from './modules/maaslin2_correlation'
include { VISUALIZATION }          from './modules/visualization'

// ── Parameter validation ───────────────────────────────────────────────────
def validate_params() {
    if (!params.samplesheet) {
        error "ERROR: --samplesheet is required (CSV with sample_id,fastq columns)"
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
    if (!params.metadata) {
        log.warn "No --metadata provided. Downstream analyses (diversity, LEfSe, MaAsLin2) will be skipped."
    }
}

// ── Log pipeline info ──────────────────────────────────────────────────────
log.info """
╔══════════════════════════════════════════════════════════════════════╗
║  WGS Metagenomics Pipeline — Nanopore Long-Read                    ║
║  ASD Biomarker Study (Multi-Sample)                                ║
╚══════════════════════════════════════════════════════════════════════╝

  Samplesheet      : ${params.samplesheet}
  Metadata         : ${params.metadata ?: 'N/A'}
  Output dir       : ${params.outdir}

  Skip host removal: ${params.skip_host_removal}
  Skip assembly    : ${params.skip_assembly}
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

    // ── Read samplesheet ────────────────────────────────────────────
    // CSV must have columns: sample_id, fastq
    ch_samples = Channel
        .fromPath(params.samplesheet, checkIfExists: true)
        .splitCsv(header: true)
        .map { row -> tuple(row.sample_id, file(row.fastq, checkIfExists: true)) }

    // ── PHASE 1: QC (per sample, parallel) ──────────────────────────
    NANOPLOT_QC(ch_samples)

    CHOPPER_FILTER(ch_samples)

    // ── PHASE 2: Host removal (per sample) ──────────────────────────
    if (!params.skip_host_removal) {
        ch_host_ref = file(params.host_reference, checkIfExists: true)
        MINIMAP2_HOST_REMOVAL(CHOPPER_FILTER.out.fastq, ch_host_ref)
        ch_clean = MINIMAP2_HOST_REMOVAL.out.fastq
    } else {
        ch_clean = CHOPPER_FILTER.out.fastq
    }

    // ── PHASE 3: Taxonomic classification (per sample) ──────────────
    ch_kraken2_db = file(params.kraken2_db, checkIfExists: true)

    KRAKEN2_CLASSIFY(ch_clean, ch_kraken2_db)

    BRACKEN_ABUNDANCE(KRAKEN2_CLASSIFY.out.report, ch_kraken2_db)

    KRONA_PLOT(KRAKEN2_CLASSIFY.out.report)

    // ── PHASE 4: Metagenome assembly + functional profiling (per sample)
    if (!params.skip_functional) {

        if (!params.skip_assembly) {
            // Assembly route (recommended): metaFlye → QUAST → filter → Prodigal
            METAFLYE_ASSEMBLE(ch_clean)
            QUAST_ASSESSMENT(METAFLYE_ASSEMBLE.out.assembly)
            SEQKIT_FILTER(METAFLYE_ASSEMBLE.out.assembly)
            PRODIGAL_PREDICT(SEQKIT_FILTER.out.fasta)
        } else {
            // Fallback: run Prodigal directly on reads (no assembly)
            FASTQ_TO_FASTA(ch_clean)
            PRODIGAL_PREDICT(FASTQ_TO_FASTA.out.fasta)
        }

        ch_eggnog_db = file(params.eggnog_db, checkIfExists: true)
        EGGNOG_ANNOTATE(PRODIGAL_PREDICT.out.proteins, ch_eggnog_db)

        EGGNOG_TO_PATHWAY_MATRIX(EGGNOG_ANNOTATE.out.annotations)

        // Collect all pathway files and merge
        ch_pathway_files = EGGNOG_TO_PATHWAY_MATRIX.out.pathway_abundance
            .map { sample_id, path -> path }
            .collect()
        MERGE_PATHWAYS(ch_pathway_files, 'pathways')
    }

    // ── PHASE 5: Merge per-sample Bracken outputs ───────────────────
    ch_bracken_files = BRACKEN_ABUNDANCE.out.abundance
        .map { sample_id, path -> path }
        .collect()
    MERGE_TAXA(ch_bracken_files, 'taxa')

    // ── PHASE 6: Diversity analysis (cohort-level) ──────────────────
    if (!params.skip_diversity && params.metadata) {
        ch_metadata = file(params.metadata, checkIfExists: true)
        DIVERSITY_ANALYSIS(MERGE_TAXA.out.matrix, ch_metadata)
    }

    // ── PHASE 7: LEfSe differential enrichment (cohort-level) ──────
    if (!params.skip_lefse && params.metadata) {
        ch_metadata_lefse = file(params.metadata, checkIfExists: true)

        LEFSE_TAXA(MERGE_TAXA.out.matrix, ch_metadata_lefse, 'taxa')

        if (!params.skip_functional) {
            LEFSE_PATHWAYS(MERGE_PATHWAYS.out.matrix, ch_metadata_lefse, 'pathways')
        }
    }

    // ── PHASE 8: MaAsLin2 biomarker-microbiome correlation (cohort) ─
    if (!params.skip_maaslin2 && params.metadata) {
        ch_metadata_maaslin = file(params.metadata, checkIfExists: true)
        MAASLIN2_ANALYSIS(MERGE_TAXA.out.matrix, ch_metadata_maaslin)
    }

    // ── PHASE 9: Visualization (cohort-level) ───────────────────────
    if (!params.skip_visualization && params.metadata) {
        ch_metadata_viz = file(params.metadata, checkIfExists: true)

        // Collect results from upstream steps (use placeholder if skipped)
        div_results = (!params.skip_diversity && params.metadata) ?
            DIVERSITY_ANALYSIS.out.results :
            file("${projectDir}/assets/PLACEHOLDER")
        maa_results = (!params.skip_maaslin2 && params.metadata) ?
            MAASLIN2_ANALYSIS.out.results :
            file("${projectDir}/assets/PLACEHOLDER")

        VISUALIZATION(
            MERGE_TAXA.out.matrix,
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
