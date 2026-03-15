process VISUALIZATION {
    label 'process_medium'
    container 'bio-r-metagenomics:1.0.0'
    publishDir "${params.outdir}/09_visualization", mode: 'copy'

    input:
    path merged_abundance
    path metadata
    path diversity_results
    path maaslin2_results

    output:
    path "figures/", emit: figures

    script:
    """
    mkdir -p figures

    Rscript ${projectDir}/bin/visualization.R \\
        --abundance ${merged_abundance} \\
        --metadata ${metadata} \\
        --diversity_dir ${diversity_results} \\
        --maaslin2_dir ${maaslin2_results} \\
        --outdir figures
    """
}
