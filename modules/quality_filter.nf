process CHOPPER_FILTER {
    tag "${sample_id}"
    label 'process_medium'
    container 'quay.io/biocontainers/chopper:0.9.0--hdfd78af_0'

    input:
    val  sample_id
    path fastq

    output:
    path "${sample_id}_filtered.fastq", emit: fastq

    script:
    """
    chopper \\
        -q ${params.min_quality} \\
        -l ${params.min_length} \\
        --threads ${task.cpus} \\
        -i ${fastq} \\
        > ${sample_id}_filtered.fastq
    """
}
