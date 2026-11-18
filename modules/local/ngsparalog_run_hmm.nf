process RUN_HMM {
    tag "segment_${segment_id}"
    label 'process_medium'

    container "${projectDir}/containers/ngsparalog.sif"
    
    beforeScript = '''
    module load singularity-ce/singularity-ce.4.1.4
    export SINGULARITY_TMPDIR=${params.singularity_tmpdir}
    export TMPDIR=${params.tmpdir}
    mkdir -p "$SINGULARITY_TMPDIR"
    '''

    input:
    tuple val(segment_id), path(lr_file), path(cov_file)
    path paramfile

    output:
    tuple val(segment_id), path("Mpyr_hmm_region_${segment_id}.rf"), emit: hmm_regions

    script:
    """
    source activate dupHMM

    # Filter coverage to just LR sites
    awk 'NR==FNR {key[\$1"\\t"\$2]=1; next} (\$1"\\t"\$2) in key' \\
        ${lr_file} \\
        ${cov_file} \\
        > average_ind_coverage_filtered_${segment_id}.txt

    Rscript /opt/ngsParalog/dupHMM.R \\
        --lrfile ${lr_file} \\
        --outfile Mpyr_hmm_region_${segment_id} \\
        --covfile average_ind_coverage_filtered_${segment_id}.txt \\
        --lrquantile ${params.lr_quantile} \\
        --paramfile=${paramfile}
    """
}