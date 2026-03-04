process MAASLIN2_ANALYSIS {
    tag "${sample_id}"
    label 'process_medium'
    container 'bio-r-metagenomics:1.0.0'
    publishDir "${params.outdir}/07_maaslin2", mode: 'copy'

    input:
    val  sample_id
    path bracken_abundance
    path metadata

    output:
    path "${sample_id}_maaslin2_results/", emit: results

    script:
    """
    mkdir -p ${sample_id}_maaslin2_results

    Rscript ${projectDir}/bin/maaslin2_analysis.R \\
        --abundance ${bracken_abundance} \\
        --metadata ${metadata} \\
        --outdir ${sample_id}_maaslin2_results \\
        --fixed_effects ${params.fixed_effects} \\
        --reference_level ${params.reference_level} \\
        --sample_id ${sample_id}
    """
}
