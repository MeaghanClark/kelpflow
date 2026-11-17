process FILTER_BY_COVERAGE {
    label 'process_low'

    input: 
    path stats_tsv

    output:
    path "samples_passing_coverage.txt", emit: pass_cov_samples
    path "samples_failing_coverage.txt", emit: fail_cov_samples
    
    script:
    """
    awk -F'\\t' -v min_cov=${params.min_cov} \\
        'NR==1 {next} \$5+0 > min_cov {print \$1 > "samples_passing_coverage.txt"; next} {print \$1 > "samples_failing_coverage.txt"}' \\
        ${stats_tsv}

    touch samples_passing_coverage.txt samples_failing_coverage.txt
    """
}