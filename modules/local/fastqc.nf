process FASTQC{
    tag "$meta.id"
    label 'process_low'

    module 'fastqc/0.11.9'

    publishDir "${params.outdir}/read_QC", mode: 'copy'

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("*.fastqc.html"), emit: html
    tuple val(meta), path("*.fastqc.zip"), emit: zip

    script:
    """
    fastqc --threads ${task.cpus}  ${reads}
    """
}
