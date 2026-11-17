process TRAIN_HMM {
    label 'process_medium'

    container "${projectDir}/containers/ngsparalog.sif"
    beforeScript = '''
    module load singularity-ce/singularity-ce.4.1.4
    export SINGULARITY_TMPDIR=/hb/scratch/miclark/singularity_tmp
    export TMPDIR=/hb/scratch/miclark/singularity_tmp
    mkdir -p "$SINGULARITY_TMPDIR"
    '''

    input:
    path lr_file
    path cov_file

    output:
    path "Mpyr_hmm_params.par", emit: paramfile

    script:
    """
    source activate dupHMM

    # Filter coverage file to LR sites
    awk 'NR==FNR {key[\$1"\\t"\$2]=1; next} (\$1"\\t"\$2) in key' \\
        ${lr_file} \\
        ${cov_file} \\
        > average_ind_coverage_filtered_all.txt


    Rscript /opt/ngsParalog/dupHMM.R \\
        --lrfile ${lr_file} \\
        --outfile Mpyr_hmm_params \\
        --covfile average_ind_coverage_filtered_all.txt \\
        --lrquantile ${params.lr_quantile} \\
        --paramOnly=1
    """
}