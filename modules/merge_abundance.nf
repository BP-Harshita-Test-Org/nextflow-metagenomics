process MERGE_ABUNDANCE {
    label 'process_low'
    container 'python:3.11-slim'
    publishDir "${params.outdir}/05_merged", mode: 'copy'

    input:
    path abundance_files
    val  analysis_type    // 'taxa' or 'pathways'

    output:
    path "merged_${analysis_type}_abundance.tsv", emit: matrix

    script:
    """
    python ${projectDir}/bin/merge_abundance.py \\
        --input_files ${abundance_files} \\
        --output merged_${analysis_type}_abundance.tsv \\
        --type ${analysis_type}
    """
}
