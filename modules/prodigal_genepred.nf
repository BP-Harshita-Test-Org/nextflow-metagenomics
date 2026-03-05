process PRODIGAL_PREDICT {
    tag "${sample_id}"
    label 'process_medium'
    container 'nanozoo/prodigal:2.6.3--2769024'
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
    # Convert FASTQ to FASTA (Prodigal requires FASTA input)
    sed -n '1~4s/^@/>/p;2~4p' ${fastq} > ${sample_id}_input.fasta

    prodigal \\
        -i ${sample_id}_input.fasta \\
        -a ${sample_id}_proteins.faa \\
        -o ${sample_id}_genes.gff \\
        -d ${sample_id}_genes.fna \\
        -p meta \\
        -f gff
    """
}
