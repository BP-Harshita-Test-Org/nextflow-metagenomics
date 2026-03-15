process FASTQ_TO_FASTA {
    tag "${sample_id}"
    label 'process_low'
    container 'nanozoo/seqtk:1.3--dc0d16b'

    input:
    tuple val(sample_id), path(fastq)

    output:
    tuple val(sample_id), path("${sample_id}_input.fasta"), emit: fasta

    script:
    """
    seqtk seq -a ${fastq} > ${sample_id}_input.fasta
    """
}

process PRODIGAL_PREDICT {
    tag "${sample_id}"
    label 'process_medium'
    container 'hyattpd/prodigal:v2.6.3'
    publishDir "${params.outdir}/04_functional/prodigal", mode: 'copy'

    input:
    tuple val(sample_id), path(fasta)

    output:
    tuple val(sample_id), path("${sample_id}_proteins.faa"), emit: proteins
    path "${sample_id}_genes.gff", emit: gff
    path "${sample_id}_genes.fna", emit: nucleotides

    script:
    """
    prodigal \\
        -i ${fasta} \\
        -a ${sample_id}_proteins.faa \\
        -o ${sample_id}_genes.gff \\
        -d ${sample_id}_genes.fna \\
        -p meta \\
        -f gff
    """
}
