process TRIM_READS{
    tag "$meta.id"
    label 'process_medium'

    module  'fastp/0.23.2'

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("*_r1.fq.gz"),      emit: forward
    tuple val(meta), path("*_r2.fq.gz"),      emit: reverse
    tuple val(meta), path("*_merged.fq.gz"),  emit: merged
    tuple val(meta), path("*_unprd_*.fq.gz"), emit: unpaired
    tuple val(meta), path("*_fail.fq.gz"),    emit: failed
    tuple val(meta), path("*.html"),          emit: html
    tuple val(meta), path("*.json"),          emit: json
    
    script:
    """
    fastp \
	    --in1 ${reads[0]} \
	    --in2 ${reads[1]} \
	    --out1 ${meta.id}_r1.fq.gz \
	    --out2 ${meta.id}_r2.fq.gz \
	    --unpaired1 ${meta.id}_unprd_1.fq.gz \
	    --unpaired2 ${meta.id}_unprd_2.fq.gz \
	    --failed_out ${meta.id}_fail.fq.gz \
	    --merge \
	    --merged_out ${meta.id}_merged.fq.gz \
	    -h ${meta.id}_fastp.html \
	    -j ${meta.id}_fastp.json \
	    --qualified_quality_phred 20 \
	    --unqualified_percent_limit 40 \
	    --length_required 30 \
	    --detect_adapter_for_pe \
	    --correction \
	    --overlap_len_require 10 \
	    --cut_tail \
	    --cut_tail_window_size 4 \
	    --cut_tail_mean_quality 15 \
	    --trim_poly_g \
	    --poly_g_min_len 10 \
	    --report_title "read_trim" 
    """
}
