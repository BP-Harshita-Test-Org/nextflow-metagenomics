process FASTQ_TO_FASTA {
    tag "${sample_id}"
    label 'process_low'
    container 'nanozoo/seqtk:1.4--48b1546'

    input:
    val  sample_id
    path fastq

    output:
    path "${sample_id}_input.fasta", emit: fasta

    script:
    """
    seqtk seq -a ${fastq} > ${sample_id}_input.fasta
    """
}

process PRODIGAL_PREDICT {
    tag "${sample_id}"
    label 'process_medium'
    container 'nanozoo/prodigal:2.6.3--2769024'
    publishDir "${params.outdir}/04_functional/prodigal", mode: 'copy'

    input:
    val  sample_id
    path fasta

    output:
    path "${sample_id}_proteins.faa", emit: proteins
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
