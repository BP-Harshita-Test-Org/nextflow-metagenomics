process VISUALIZATION {
    tag "${sample_id}"
    label 'process_medium'
    container 'bio-r-metagenomics:1.0.0'
    publishDir "${params.outdir}/08_visualization", mode: 'copy'

    input:
    val  sample_id
    path bracken_abundance
    path metadata
    path diversity_results
    path maaslin2_results

    output:
    path "${sample_id}_figures/", emit: figures

    script:
    """
    mkdir -p ${sample_id}_figures

    Rscript ${projectDir}/bin/visualization.R \\
        --abundance ${bracken_abundance} \\
        --metadata ${metadata} \\
        --diversity_dir ${diversity_results} \\
        --maaslin2_dir ${maaslin2_results} \\
        --outdir ${sample_id}_figures \\
        --sample_id ${sample_id}
    """
}
