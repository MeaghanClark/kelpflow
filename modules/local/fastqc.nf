process FASTQC{
    tag "$meta.id"
    label 'process_low'

    module 'fastqc/0.12.1'


    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("*_fastqc.html"), emit: html
    tuple val(meta), path("*_fastqc.zip"), emit: zip

    script:
    """
    fastqc --threads ${task.cpus} ${reads}
    """
}
