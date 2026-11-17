process RUN_PCANGSD {
    label 'process_high'
    module 'miniconda3'

    input:
    path(beagle)

    output:
    path("pca.cov"),  emit: cov
    path("pca.*"),    emit: all

    script:
    """
    set -eo pipefail 
    
    eval "\$(conda shell.bash hook)"
    conda activate pcangsd
    
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