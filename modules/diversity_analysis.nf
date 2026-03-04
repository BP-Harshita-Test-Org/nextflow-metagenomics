process DIVERSITY_ANALYSIS {
    tag "${sample_id}"
    label 'process_medium'
    container 'bio-r-metagenomics:1.0.0'
    publishDir "${params.outdir}/05_diversity", mode: 'copy'

    input:
    val  sample_id
    path bracken_abundance
    path metadata

    output:
    path "${sample_id}_diversity_results/", emit: results

    script:
    """
    mkdir -p ${sample_id}_diversity_results

    Rscript ${projectDir}/bin/diversity_analysis.R \\
        --abundance ${bracken_abundance} \\
        --metadata ${metadata} \\
        --outdir ${sample_id}_diversity_results \\
        --fdr_alpha ${params.fdr_alpha} \\
        --fdr_ancombc ${params.fdr_ancombc} \\
        --sample_id ${sample_id}
    """
}
