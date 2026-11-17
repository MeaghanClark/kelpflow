process PARSE_QUALIMAP {
    tag "$meta.id"
    label 'process_low'

    input:
    tuple val(meta), path(qualimap_dir)

    output:
    tuple val(meta), path("${meta.id}_qualimap_stats.tsv"), emit: tsv

    script:
    """
    results_file="${qualimap_dir}/genome_results.txt"

    # desired stats: total reads, mapped reads, percent mapped, mean coverage, sd of coverage, mean mapping quality, error rate
    total_reads=\$(grep "number of reads" "\$results_file" | grep -o '[0-9,]*' | tr -d ',')
    mapped_reads=\$(grep "number of mapped reads" "\$results_file" | grep -oP '\\d[\\d,]*(?= \\()' | tr -d ',')
    pct_mapped=\$(grep "number of mapped reads" "\$results_file" | grep -oP '(?<=\\()[0-9.]+(?=%\\))')
    mean_coverage=\$(grep "mean coverageData" "\$results_file" | grep -oP '[0-9.]+(?=X)')
    std_coverage=\$(grep "std coverageData" "\$results_file" | grep -oP '[0-9.]+(?=X)')
    mean_mapq=\$(grep "mean mapping quality" "\$results_file" | grep -oP '[0-9.]+')
    error_rate=\$(grep "error rate" "\$results_file" | grep -oP '[0-9.]+')
    gc_pct=\$(grep "GC percentage" "\$results_file" | grep -oP '[0-9.]+(?=%)')

    echo -e "sample\\ttotal_reads\\tmapped_reads\\tpct_mapped\\tmean_coverage\\tstd_coverage\\tmean_mapq\\terror_rate\\tgc_pct" > ${meta.id}_qualimap_stats.tsv
    echo -e "${meta.id}\\t\${total_reads}\\t\${mapped_reads}\\t\${pct_mapped}\\t\${mean_coverage}\\t\${std_coverage}\\t\${mean_mapq}\\t\${error_rate}\\t\${gc_pct}" >> ${meta.id}_qualimap_stats.tsv
    """
}