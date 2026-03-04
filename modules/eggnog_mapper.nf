process EGGNOG_ANNOTATE {
    tag "${sample_id}"
    label 'process_high'
    container 'nanozoo/eggnog-mapper:2.1.12--4f2b6c0'
    publishDir "${params.outdir}/04_functional/eggnog", mode: 'copy'

    input:
    val  sample_id
    path proteins
    path eggnog_db

    output:
    path "${sample_id}_eggnog.emapper.annotations", emit: annotations
    path "${sample_id}_eggnog.emapper.seed_orthologs", emit: orthologs

    script:
    """
    emapper.py \\
        -i ${proteins} \\
        --output ${sample_id}_eggnog \\
        --data_dir ${eggnog_db} \\
        -m diamond \\
        --cpu ${task.cpus} \\
        --override \\
        --go_evidence all
    """
}
