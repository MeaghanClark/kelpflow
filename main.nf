nextflow.enable.dsl = 2

// Import modules
include { FASTQC        } from './modules/local/fastqc.nf'
include { MULTIQC       } from './modules/local/multiqc.nf'

// Workflow

workflow {
    Channel
        .fromPath("${params.input_dir}/*.fastq.gz")
        .map { f -> tuple([id: f.simpleName], f) }
        .set { ch_reads }

    FASTQC(ch_reads)

    ch_multiqc_files = FASTQC.out.zip
        .map { meta, zip -> zip }
        .collect()

    MULTIQC(ch_multiqc_files, params.report_name)
}
