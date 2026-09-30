process MULTIQC {
    label 'process_medium'

    module 'multiqc/1.27'

    input:
    path 'fastqc/*'
    val  report_name

    output:
    path "*.html",                emit: report
    path "${report_name}_data",   emit: data
    path "*_plots",               emit: plots, optional: true

    script:
    """
    multiqc \\
        -v -p -ip -f \\
        --data-dir \\
        --data-format tsv \\
        --cl-config "max_table_rows: 3000" \\
        --filename ${report_name} \\
        .
    """
}