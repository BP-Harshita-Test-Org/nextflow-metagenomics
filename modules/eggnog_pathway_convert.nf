process EGGNOG_TO_PATHWAY_MATRIX {
    tag "${sample_id}"
    label 'process_low'
    // Uses the same R metagenomics image which has Python available,
    // or swap for any lightweight Python image
    container 'python:3.11-slim'

    input:
    val  sample_id
    path annotations   // *.emapper.annotations from EGGNOG_ANNOTATE

    output:
    path "${sample_id}_pathway_abundance.tsv", emit: pathway_abundance

    script:
    """
    python ${projectDir}/bin/eggnog_to_pathway_matrix.py \\
        --annotations ${annotations} \\
        --output ${sample_id}_pathway_abundance.tsv \\
        --sample_id ${sample_id}
    """
}
