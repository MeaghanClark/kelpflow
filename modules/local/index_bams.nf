process INDEX_BAMS {
    tag "$meta.id"
    label 'process_high'

    module 'samtools/1.21'

    input: 
    tuple val(meta), path(bam)*

    output: 
    tuple val(meta), path("*.sorted.bam"), path("*.sorted.bam.bai"), emit: bam_bai

    script:
    """
    samtools sort -@ ${task.cpus} ${bam} -o ${meta.id}.sorted.bam
    samtools index -@ ${task.cpus} ${meta.id}.sorted.bam
    """
}
