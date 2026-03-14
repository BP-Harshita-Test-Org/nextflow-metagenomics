process BRACKEN_ABUNDANCE {
    tag "${sample_id}"
    label 'process_medium'
    container 'staphb/bracken:latest'
    publishDir "${params.outdir}/03_taxonomy/bracken", mode: 'copy'

    input:
    val  sample_id
    path kraken2_report
    path kraken2_db

    output:
    path "${sample_id}_bracken_species.tsv", emit: abundance
    path "${sample_id}_bracken_report.txt", emit: report

    script:
    """
    bracken \\
        -d ${kraken2_db} \\
        -i ${kraken2_report} \\
        -o ${sample_id}_bracken_species.tsv \\
        -w ${sample_id}_bracken_report.txt \\
        -r ${params.bracken_read_length} \\
        -l ${params.bracken_level} \\
        -t ${params.bracken_threshold}
    """
}
