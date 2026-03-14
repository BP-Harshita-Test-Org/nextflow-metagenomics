process LEFSE_ANALYSIS {
    tag "${sample_id}"
    label 'process_medium'
    container 'biobakery/lefse:latest'
    publishDir "${params.outdir}/06_lefse", mode: 'copy'

    input:
    val  sample_id
    path abundance_table
    path metadata
    val  analysis_type   // 'taxa' or 'pathways'

    output:
    path "${sample_id}_lefse_${analysis_type}_results.tsv", emit: results
    path "${sample_id}_lefse_${analysis_type}_plot.pdf", optional: true, emit: plot

    script:
    """
    # Format input for LEfSe
    python ${projectDir}/bin/format_lefse_input.py \\
        --abundance ${abundance_table} \\
        --metadata ${metadata} \\
        --output lefse_input.txt

    # Run LEfSe
    lefse_format_input.py lefse_input.txt lefse_formatted.in -c 1 -u 2 -o 1000000
    lefse_run.py lefse_formatted.in ${sample_id}_lefse_${analysis_type}_results.tsv -l ${params.lda_threshold}
    lefse_plot_res.py ${sample_id}_lefse_${analysis_type}_results.tsv ${sample_id}_lefse_${analysis_type}_plot.pdf --format pdf || true
    """
}
