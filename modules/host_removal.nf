process MINIMAP2_HOST_REMOVAL {
    tag "${sample_id}"
    label 'process_high'
    container 'nanozoo/minimap2:2.28--9e3bd01'
    publishDir "${params.outdir}/02_host_removal", mode: 'copy', pattern: '*_host_stats.txt'

    input:
    val  sample_id
    path fastq
    path host_ref

    output:
    path "${sample_id}_nonhost.fastq", emit: fastq
    path "${sample_id}_host_stats.txt", emit: stats

    script:
    """
    minimap2 \\
        -a -x map-ont \\
        -t ${task.cpus} \\
        ${host_ref} \\
        ${fastq} \\
        | samtools view -bS -@ ${task.cpus} - \\
        | samtools sort -@ ${task.cpus} -o aligned.bam

    samtools index aligned.bam

    # Extract unmapped (non-host) reads
    samtools view -b -f 4 aligned.bam \\
        | samtools fastq - > ${sample_id}_nonhost.fastq

    # Host removal stats
    TOTAL=\$(samtools view -c aligned.bam)
    HOST=\$(samtools view -c -F 4 aligned.bam)
    NONHOST=\$(samtools view -c -f 4 aligned.bam)
    echo "total_reads: \${TOTAL}" > ${sample_id}_host_stats.txt
    echo "host_reads: \${HOST}" >> ${sample_id}_host_stats.txt
    echo "nonhost_reads: \${NONHOST}" >> ${sample_id}_host_stats.txt
    echo "host_fraction: \$(awk "BEGIN {printf \\"%.4f\\", \${HOST}/\${TOTAL}}")" >> ${sample_id}_host_stats.txt
    """
}
