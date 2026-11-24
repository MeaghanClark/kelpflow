process GENERATE_PCA_LIKELIHOODS {
    label 'process_high'
    module 'angsd/0.940'

    input:
    path(bamlist)
    path(bam_files)
    path(bai_files)
    path(sites)          // pruned sites from PRUNE_LD_GRAPH
    path(depth_cutoff)   // depth cutoff from DEPTH_HIST

    output:
    path("pca_input.beagle.gz"), emit: beagle
    path("pca_input.mafs.gz"),   emit: mafs

    script:
    def max_depth_val = "cat ${depth_cutoff} | awk '{print \$NF}'"

    """
    set -euo pipefail

    # ANGSD -sites wants chr<TAB>pos. Convert chr:pos if the pruned list still has colons.
    if grep -q ':' ${sites}; then
      sed 's/:\\([^:]*\\)\$/\\t\\1/' ${sites} > sites.tsv
    else
      cp ${sites} sites.tsv
    fi
    angsd sites index sites.tsv

    NINDS=\$(wc -l < ${bamlist})
    MIN_IND=\$(awk -v n=\$NINDS -v f=${params.pca_min_ind_frac} 'BEGIN{printf "%d", n*f+0.5}')

    angsd \\
      -b ${bamlist} \\
      -sites sites.tsv \\
      -GL 1 \\
      -doGlf 2 \\
      -doMajorMinor 1 \\
      -doMaf 1 \\
      -doCounts 1 \\
      -SNP_pval ${params.pca_snp_pval} \\
      -minInd \$MIN_IND \\
      -minMapQ 20 -minQ 20 \\
      -setmaxDepth ${max_depth_val} \\
      -nThreads ${task.cpus} \\
      -out pca_input

    N=\$(zcat pca_input.beagle.gz | tail -n +2 | wc -l)
    [ "\$N" -gt 0 ] || { echo "ERROR: 0 sites in beagle" >&2; exit 1; }
    echo "sites for PCA: \$N"
    """
}