process KRAKEN2_CLASSIFY {
    tag "${sample_id}"
    label 'process_high_memory'
    container 'staphb/kraken2:latest'
    publishDir "${params.outdir}/03_taxonomy/kraken2", mode: 'copy'

    input:
    tuple val(sample_id), path(fastq)
    path kraken2_db

    output:
    path "${sample_id}_kraken2_output.txt", emit: output
    tuple val(sample_id), path("${sample_id}_kraken2_report.txt"), emit: report

    script:
    """
    kraken2 \\
        --db ${kraken2_db} \\
        --threads ${task.cpus} \\
        --output ${sample_id}_kraken2_output.txt \\
        --report ${sample_id}_kraken2_report.txt \\
        --use-names \\
        --memory-mapping \\
        ${fastq}
    """
}
