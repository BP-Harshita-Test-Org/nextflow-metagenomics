process METAFLYE_ASSEMBLE {
    tag "${sample_id}"
    label 'process_high'
    container 'nanozoo/flye:2.9.3'
    publishDir "${params.outdir}/04_functional/assembly", mode: 'copy'

    input:
    tuple val(sample_id), path(fastq)

    output:
    tuple val(sample_id), path("${sample_id}_assembly.fasta"), emit: assembly
    path "${sample_id}_assembly_info.txt", emit: info

    script:
    """
    flye \\
        --nano-raw ${fastq} \\
        --out-dir assembly/ \\
        --meta \\
        --threads ${task.cpus}

    cp assembly/assembly.fasta ${sample_id}_assembly.fasta
    cp assembly/assembly_info.txt ${sample_id}_assembly_info.txt
    """
}

process QUAST_ASSESSMENT {
    tag "${sample_id}"
    label 'process_low'
    container 'staphb/quast:latest'
    publishDir "${params.outdir}/04_functional/quast", mode: 'copy'

    input:
    tuple val(sample_id), path(assembly)

    output:
    path "${sample_id}_quast_results/", emit: results

    script:
    """
    quast ${assembly} \\
        -o ${sample_id}_quast_results \\
        --min-contig ${params.min_contig_length} \\
        --threads ${task.cpus}
    """
}

process SEQKIT_FILTER {
    tag "${sample_id}"
    label 'process_low'
    container 'quay.io/biocontainers/seqkit:2.6.1--h9ee0642_0'

    input:
    tuple val(sample_id), path(assembly)

    output:
    tuple val(sample_id), path("${sample_id}_filtered_contigs.fasta"), emit: fasta

    script:
    """
    seqkit seq -m ${params.min_contig_length} ${assembly} > ${sample_id}_filtered_contigs.fasta
    """
}
