process FILTER_PARALOGS {
    module 'bedtools/2.26.0'

    input:
    path(putative_sites)
    path(hmm_regions_files)

    output:
    path("filtered_sites.txt"), emit: sites

    script:
    """
    set -euo pipefail
    
    cat ${hmm_regions_files} | sort -k1,1 -k2,2n > problem_regions.bed
    
    # collectFile order is arbitrary, so sort; dedupe in case segments overlap
    sort -k1,1 -k2,2n -u ${putative_sites} \\
      | awk 'BEGIN{OFS="\\t"} {print \$1, \$2-1, \$2}' > all_sites.bed
    
    bedtools intersect -a all_sites.bed -b problem_regions.bed -v | cut -f1,3 > filtered_sites.txt
    
    [ -s filtered_sites.txt ] || { echo "ERROR: 0 sites after paralog filter" >&2; exit 1; }
    """
}