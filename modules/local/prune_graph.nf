process PRUNE_LD_GRAPH {
    label 'process_high'

    input:
    tuple val(meta), path(ld_file)

    output:
    tuple val(meta), path("pruned_safe_sites.txt"), emit: pruned_sites
    path("ld_pruning_report.txt"), emit: report

    script:
    """
    set -euo pipefail
    
    PRUNEGRAPH="${params.prune_graph_path}"
    
    \$PRUNEGRAPH \\
      --in ${ld_file} \\
      --weight "column_7" \\
      --filter "column_3 <= ${params.ld_max_bp_dist} && column_7 >= ${params.ld_min_r2}" \\
      --n-threads ${task.cpus} \\
      --out unlinked_sites.pos
    
    # chr:pos -> chr<TAB>pos, sorted (split on the LAST colon in case contig names contain one)
    awk 'BEGIN{OFS="\\t"} {i=match(\$0,/:[^:]*\$/); print substr(\$0,1,i-1), substr(\$0,i+1)}' unlinked_sites.pos \\
      | sort -k1,1 -k2,2n > pruned_safe_sites.txt
    
    [ -s pruned_safe_sites.txt ] || { echo "ERROR: pruning produced 0 sites" >&2; exit 1; }
    
    {
      echo "LD Pruning Report"
      echo "================"
      echo "Input LD file: ${ld_file}"
      echo "Max distance: ${params.ld_max_bp_dist} bp"
      echo "Min r2: ${params.ld_min_r2}"
      echo "Output sites: \$(wc -l < pruned_safe_sites.txt)"
    } > ld_pruning_report.txt
    """
}