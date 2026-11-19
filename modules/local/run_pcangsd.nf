process RUN_PCANGSD {
    label 'process_high'

    module 'miniconda3'
    conda 'bioconda::pcangsd=1.36.4'

    input:
    path(beagle)

    output:
    path("pca.cov"),  emit: cov
    path("pca.*"),    emit: all

    script:
    """
    set -eo pipefail 
        
    pcangsd \\
      --beagle ${beagle} \\
      --maf ${params.pca_maf} \\
      --pcadapt \\
      --inbreed-samples \\
      --threads ${task.cpus} \\
      --sites-save \\
      --out pca
    """
}