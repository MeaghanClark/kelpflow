process MAKE_NGSPARALOG_INPUT {
    tag "segment_${segment_id}"
    label 'process_medium'

    // module 'angsd/0.940:samtools/1.21'

    input:
    tuple val(segment_id), path(segment_file)
    path bamlist
    path bam_files
    path bai_files

    output:
    tuple val(segment_id), path("putative_varsites_${segment_id}.pos"),          emit: varsites
    tuple val(segment_id), path("average_ind_coverage_${segment_id}.txt"),       emit: coverage
    tuple val(segment_id), path("total_ind_coverage_${segment_id}.txt"),         emit: coverage_total
    tuple val(segment_id), path("putative_snps_${segment_id}.mafs.gz"),          emit: mafs

    script:
    """
    source /opt/ohpc/admin/lmod/lmod/init/bash
    
    NINDS=\$(wc -l < ${bamlist})
    MIN_IND=\$(( (NINDS + 1) / 2 ))

    module load angsd/0.940
    angsd \\
        -b ${bamlist} \\
        -GL 1 \\
        -doMaf 1 \\
        -doMajorMinor 1 \\
        -SNP_pval 1e-5 \\
        -minInd \$MIN_IND \\
        -nThreads ${task.cpus} \\
        -rf ${segment_file} \\
        -out putative_snps_${segment_id}

    zcat putative_snps_${segment_id}.mafs.gz | cut -f 1,2 | tail -n +2 \\
        > putative_varsites_${segment_id}.pos

    module load samtools/1.21
    samtools depth -a -f ${bamlist} -b putative_varsites_${segment_id}.pos \\
        | awk -v avg_out="average_ind_coverage_${segment_id}.txt" \\
              -v tot_out="total_ind_coverage_${segment_id}.txt" \\
              'BEGIN{OFS="\\t"} {
                  sum=0
                  for(i=3; i<=NF; i++) sum+=\$i
                  print \$1, \$2, sum/(NF-2) > avg_out
                  print \$1, \$2, sum        > tot_out
              }'

    """
}