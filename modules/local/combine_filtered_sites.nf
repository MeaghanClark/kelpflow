process COMBINE_FILTERED_SITES {
    input:
    tuple val(meta), path(pruned_sites)
    tuple val(meta2), path(safe_sites)  // sites not in problem regions

    output:
    tuple val(meta), path("pruned_safe_sites.txt"), emit: final_sites

    script:
    """
    # Sites that are: (1) not in problem regions AND (2) not in high LD
    # Get intersection (sites must appear in both files)
    sort ${pruned_sites} > pruned_sorted.txt
    zcat ${safe_sites} | sort > safe_sorted.txt

    comm -12 pruned_sorted.txt safe_sorted.txt > pruned_safe_sites.txt

    echo "Final site count: \$(wc -l < pruned_safe_sites.txt)"
    """
}