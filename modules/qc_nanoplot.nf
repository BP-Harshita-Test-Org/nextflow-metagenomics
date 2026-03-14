process NANOPLOT_QC {
    tag "${sample_id}"
    label 'process_low'
    container 'staphb/nanoplot:latest'
    publishDir "${params.outdir}/01_qc/nanoplot", mode: 'copy'

    input:
    val  sample_id
    path fastq

    output:
    path "${sample_id}_nanoplot/", emit: report_dir
    path "${sample_id}_nanoplot/NanoStats.txt", emit: stats

    script:
    """
    NanoPlot \\
        --fastq ${fastq} \\
        -o ${sample_id}_nanoplot \\
        -t ${task.cpus} \\
        --plots dot
    """
}
