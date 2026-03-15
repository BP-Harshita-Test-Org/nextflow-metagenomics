process EGGNOG_TO_PATHWAY_MATRIX {
    tag "${sample_id}"
    label 'process_low'
    container 'python:3.11-slim'

    input:
    tuple val(sample_id), path(annotations)

    output:
    tuple val(sample_id), path("${sample_id}_pathway_abundance.tsv"), emit: pathway_abundance

    script:
    """
    python ${projectDir}/bin/eggnog_to_pathway_matrix.py \\
        --annotations ${annotations} \\
        --output ${sample_id}_pathway_abundance.tsv \\
        --sample_id ${sample_id}
    """
}
