process MAASLIN2_ANALYSIS {
    label 'process_medium'
    container 'bio-r-metagenomics:1.0.0'
    publishDir "${params.outdir}/08_maaslin2", mode: 'copy'

    input:
    path merged_abundance
    path metadata

    output:
    path "maaslin2_results/", emit: results

    script:
    """
    mkdir -p maaslin2_results

    Rscript ${projectDir}/bin/maaslin2_analysis.R \\
        --abundance ${merged_abundance} \\
        --metadata ${metadata} \\
        --outdir maaslin2_results \\
        --fixed_effects ${params.fixed_effects} \\
        --reference_level ${params.reference_level}
    """
}
