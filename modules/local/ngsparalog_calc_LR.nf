process CALC_LR {
    tag "segment_${segment_id}"
    label 'process_medium'
    
    // NOTE: beforeScript/container are actually set via the withName selector in
    // nextflow.config (config-level directives override these), so they're omitted here
    // to avoid the two definitions silently diverging.

    input:
    tuple val(segment_id), path(region_file), path(posfile)
    path bamlist
    path bam_files
    path bai_files
    tuple path(ref_fasta), path(ref_fai)

    output:
    tuple val(segment_id), path("Mpyr_LR_${segment_id}.tsv"), emit: lr

    script:
    """
    NINDS=\$(wc -l < ${bamlist})
    MIN_IND=\$(( (NINDS + 1) / 2 ))

    samtools mpileup -b ${bamlist} -f ${ref_fasta} -R ${region_file} -l ${posfile} \\
        -A -a -q 0 -Q 0 --ff UNMAP,DUP,QCFAIL -d 110000 \\
        | ngsParalog calcLR -allow_overwrite 1 -infile - -outfile Mpyr_LR_${segment_id}.tsv \\
        -minQ 20 -minind \$MIN_IND -mincov 1
    """
}