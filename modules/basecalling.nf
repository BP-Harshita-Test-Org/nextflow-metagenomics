process DORADO_BASECALL {
    tag "${sample_id}"
    label 'process_high'
    container 'ontresearch/dorado:latest'

    input:
    val  sample_id
    path pod5_dir

    output:
    path "${sample_id}_basecalled.fastq", emit: fastq

    script:
    """
    dorado basecaller \\
        ${params.dorado_model} \\
        ${pod5_dir} \\
        > ${sample_id}_basecalled.fastq
    """
}
