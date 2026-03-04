process PRODIGAL_PREDICT {
    tag "${sample_id}"
    label 'process_medium'
    container 'biocontainers/prodigal:v2.6.3-1-deb_cv1'
    publishDir "${params.outdir}/04_functional/prodigal", mode: 'copy'

    input:
    val  sample_id
    path fastq

    output:
    path "${sample_id}_proteins.faa", emit: proteins
    path "${sample_id}_genes.gff", emit: gff
    path "${sample_id}_genes.fna", emit: nucleotides

    script:
    """
    prodigal \\
        -i ${fastq} \\
        -a ${sample_id}_proteins.faa \\
        -o ${sample_id}_genes.gff \\
        -d ${sample_id}_genes.fna \\
        -p meta \\
        -f gff
    """
}
