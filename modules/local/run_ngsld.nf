process RUN_NGSLD {
    label 'process_high'
    module 'ngsld/1.2.0'

    input:
    tuple val(meta), path(beagle), path(pos)
    
    output:
    tuple val(meta), path("*.ld"), emit: ld
    
    script:
    """    
    N_SITES=\$(zcat ${pos} | wc -l)
    N_IND=\$(( \$(zcat ${beagle} | head -1 | awk '{print NF}') / 3 - 1 ))
    
    ngsLD \\
      --geno ${beagle} \\
      --pos ${pos} \\
      --probs \\
      --n_ind \$N_IND \\
      --n_sites \$N_SITES \\
      --max_kb_dist ${params.ld_max_kb_dist} \\
      --n_threads ${task.cpus} \\
      --out genotypes_ld.ld
    """
}