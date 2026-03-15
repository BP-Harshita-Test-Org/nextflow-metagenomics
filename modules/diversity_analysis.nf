process DIVERSITY_ANALYSIS {
    label 'process_medium'
    container 'bio-r-metagenomics:1.0.0'
    publishDir "${params.outdir}/06_diversity", mode: 'copy'

    input:
    path merged_abundance
    path metadata

    output:
    path "diversity_results/", emit: results

    script:
    """
    mkdir -p diversity_results

    Rscript ${projectDir}/bin/diversity_analysis.R \\
        --abundance ${merged_abundance} \\
        --metadata ${metadata} \\
        --outdir diversity_results \\
        --fdr_alpha ${params.fdr_alpha} \\
        --fdr_ancombc ${params.fdr_ancombc}
    """
}
