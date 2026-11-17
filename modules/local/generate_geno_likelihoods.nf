process GENERATE_GENO_LIKELIHOODS {
    label 'process_high'
    module 'angsd/0.940'

    input:
    path(bamlist)
    path(bam_files)
    path(bai_files)
    path(sites)

    output:
    tuple val(meta), path("genotypes.beagle.gz"), emit: beagle
    tuple val(meta), path("ld_input.beagle.gz"), path("ld_sites.pos.gz"), emit: ld_input

    script:
    meta = [id: 'genome_wide']
    
    """
    set -euo pipefail

    angsd sites index ${sites}

    angsd \\
      -b ${bamlist} \\
      -sites ${sites} \\
      -GL 1 \\
      -doMajorMinor 1 \\
      -doMaf 1 \\
      -SNP_pval ${params.angsd_snp_pval} \\
      -minMapQ 20 -minQ 20 \\
      -doGlf 2 \\
      -nThreads ${task.cpus} \\
      -out genotypes

    # pos file from the beagle markers (contig_pos; split on the LAST underscore)
    zcat genotypes.beagle.gz | tail -n +2 | cut -f1 | \\
      awk 'BEGIN{OFS="\\t"} {n=split(\$1,a,"_"); p=a[n]; c=substr(\$1,1,length(\$1)-length(p)-1); print c,p}' | \\
      gzip > ld_sites.pos.gz

    # ngsLD geno: no header, no marker/allele columns
    zcat genotypes.beagle.gz | tail -n +2 | cut -f4- | gzip > ld_input.beagle.gz

    N_POS=\$(zcat ld_sites.pos.gz | wc -l)
    N_GENO=\$(zcat ld_input.beagle.gz | wc -l)
    [ "\$N_POS" -gt 0 ] && [ "\$N_POS" -eq "\$N_GENO" ] || { echo "ERROR: pos/geno mismatch (\$N_POS vs \$N_GENO)" >&2; exit 1; }
    """
}