process LEFSE_ANALYSIS {
    label 'process_medium'
    container 'biobakery/lefse:latest'
    publishDir "${params.outdir}/07_lefse", mode: 'copy'

    input:
    path abundance_table
    path metadata
    val  analysis_type   // 'taxa' or 'pathways'

    output:
    path "lefse_${analysis_type}_results.tsv", emit: results
    path "lefse_${analysis_type}_plot.pdf", optional: true, emit: plot

    script:
    """
    # Format input for LEfSe
    python ${projectDir}/bin/format_lefse_input.py \\
        --abundance ${abundance_table} \\
        --metadata ${metadata} \\
        --output lefse_input.txt

    # Run LEfSe
    lefse_format_input.py lefse_input.txt lefse_formatted.in -c 1 -s 2 -u 3 -o 1000000
    lefse_run.py lefse_formatted.in lefse_${analysis_type}_results.tsv -l ${params.lda_threshold}
    lefse_plot_res.py lefse_${analysis_type}_results.tsv lefse_${analysis_type}_plot.pdf --format pdf || true
    """
}
