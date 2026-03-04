process KRONA_PLOT {
    tag "${sample_id}"
    label 'process_low'
    container 'nanozoo/krona:2.7.1--e7615f7'
    publishDir "${params.outdir}/03_taxonomy/krona", mode: 'copy'

    input:
    val  sample_id
    path kraken2_report

    output:
    path "${sample_id}_krona.html", emit: html

    script:
    """
    /opt/conda/opt/krona/updateTaxonomy.sh
    /opt/conda/opt/krona/scripts/ImportTaxonomy.pl \\
        -o ${sample_id}_krona.html \\
        -t 5 -m 3 \\
        -i ${kraken2_report}
    """
}
