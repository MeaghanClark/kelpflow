process QUALIMAP{
    tag "$meta.id"
    label 'process_medium'

    module 'miniconda3/3.13'
    conda 'bioconda::qualimap'

    input:
    tuple val(meta), path(bam)
    
    output:
    tuple val(meta), path("${meta.id}_qualimap"), emit: results
    
    script:
    def javaMem = (task.memory.toGiga() * 0.75).toInteger()

    """
    qualimap bamqc \\
        -bam "${bam}" \\
        -outdir "${meta.id}_qualimap" \\
        -outformat HTML \\
        --java-mem-size=${javaMem}G \\
        -nt ${task.cpus} \\
        --cov-hist-lim 2000 \\
        --collect-overlap-pairs \\
        --paint-chromosome-limits

    echo "Done with ${meta.id}"
    """
}
