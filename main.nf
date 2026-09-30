nextflow.enable.dsl = 2

// Import modules
include { FASTQC                      } from './modules/local/fastqc.nf'
include { MULTIQC as MULTIQC_RAW      } from './modules/local/multiqc.nf'
include { MULTIQC as MULTIQC_TRIM     } from './modules/local/multiqc.nf'
include { TRIM_READS as TRIM          } from './modules/local/trim_reads.nf'
include { ALIGN                       } from './modules/local/align_to_genome.nf'
include { QUALIMAP                    } from './modules/local/qualimap.nf'
include { MULTIQC as MULTIQC_QUALIMAP } from './modules/local/multiqc.nf'
include { PARSE_QUALIMAP              } from './modules/local/parse_qualimap.nf'

// Workflow

workflow {

    // Define main workflow
    main:
    // Input channels
    // build channel for all reads in the input directory
    ch_reads = Channel.fromPath("${params.input_dir}/*.fastq.gz")
        .map { f -> tuple([id: f.simpleName], f) }
    // build paired-read channel 
    ch_read_pairs = Channel
        .fromFilePairs("${params.input_dir}/*_{1,2}.fastq.gz")
        .map { sample_id, reads -> tuple([id: sample_id], reads) }
    // build reference and associated index files channel 
    ch_ref_genome = Channel.fromPath("${params.reference_genome}*")
        .collect()
        .map{files -> 
            def fasta = files.find {it.name.endsWith('.fasta') || it.name.endsWith('.fa') } // Find the main FASTA file in the directory
            def index_only = files.findAll { !it.name.endsWith('.fasta') && !it.name.endsWith('.fa') } // Filter the array to exclude the main fasta file
            return tuple([id: fasta.baseName], fasta, index_only) // return tuple
        }

    // Start pipeline

    // Run fastQC on raw reads 
    FASTQC(ch_reads)

    // CHANNEL: Gather all FastQC zips into a single list for one MultiQC run
    ch_multiqc_files = FASTQC.out.zip
        .map { meta, zip -> zip }
        .collect()
    
    // Run MultiQC on fastQC outputs
    MULTIQC_RAW(ch_multiqc_files, params.raw_report_name)

    // Trim reads
    TRIM(ch_read_pairs)
    
    // CHANNEL: Gather all trimmed reads into a single list for one MultiQC run
    ch_trimmed_multiqc_files = TRIM.out.json // change to match trim_reads.nf output
        .map { meta, json -> json }
        .collect()

    // Run MultiQC on trimmed reads outputs
    MULTIQC_TRIM(ch_trimmed_multiqc_files, params.trimmed_report_name)

    // CHANNEL: Gather all trimmed reads into a single list for alignment
    ch_trimmed_reads = TRIM.out.forward
        .join(TRIM.out.reverse)
        .join(TRIM.out.merged)

    // CHANNEL: Alignment input
    ch_align_inputs = ch_trimmed_reads.combine(ch_ref_genome) // Combine every sample in ch_trimmed_reads with the single ch_ref_genome tuple

    // Align trimmed reads to reference genome
    ALIGN(ch_align_inputs)

    // CHANNEL: Gather filtered alignment files
    ch_filtered_bams = ALIGN.out.bam_filt

    // Run Qualimap on filtered BAM files
    QUALIMAP(ch_filtered_bams)

    // CHANNEL: Gather qualimap output for one MultiQC run
    ch_qualimap_files = QUALIMAP.out.results
        .map { meta, file -> file }
        .collect()

    // Run MultiQC on qualimap outputs
    MULTIQC_QUALIMAP(ch_qualimap_files, params.qualimap_report_name)

    // Parse qualimap statistics 
    PARSE_QUALIMAP(QUALIMAP.out.results)

    // Collect qualimap statistics and output report 
    ch_qualimap_stats = PARSE_QUALIMAP.out.tsv
        .map { meta, tsv -> tsv }
        .collectFile(
            name: 'alignment_stats.tsv', 
            keepHeader: true,
            storeDir: "${params.outdir}/summary_reports"
        )

    // Define output publishing sub-directories
        // specific path defined in output chunk beloew
    publish:
    // Read quality
        fastqc_html         = FASTQC.out.html
        fastqc_zip          = FASTQC.out.zip
        multiqc_raw_report  = MULTIQC_RAW.out.report
        multiqc_raw_data    = MULTIQC_RAW.out.data
        trim_forward        = TRIM.out.forward
        trim_reverse        = TRIM.out.reverse
        trim_merged         = TRIM.out.merged
        trim_unpaired       = TRIM.out.unpaired
        trim_failed         = TRIM.out.failed
        trim_html           = TRIM.out.html
        trim_json           = TRIM.out.json
        multiqc_trim_report = MULTIQC_TRIM.out.report
        multiqc_trim_data   = MULTIQC_TRIM.out.data
    // Alignment and QC
        bam_minfilt         = ALIGN.out.bam_minfilt
        bam_filt            = ALIGN.out.bam_filt
        qualimap_results  = QUALIMAP.out.results
        multiqc_qualimap_report = MULTIQC_QUALIMAP.out.report
        multiqc_qualimap_data   = MULTIQC_QUALIMAP.out.data
}

// Define overall publishing directory for the workflow
output {
    // Read quality
    fastqc_html              { path 'raw_read_QC/fastqc' }
    fastqc_zip               { path 'raw_read_QC/fastqc' }
    multiqc_raw_report       { path 'summary_reports' }
    multiqc_raw_data         { path 'raw_read_QC/multiqc' }
    trim_forward             { path 'trimmed_reads/paired' }
    trim_reverse             { path 'trimmed_reads/paired' }
    trim_merged              { path 'trimmed_reads/merged' }
    trim_unpaired            { path 'trimmed_reads/unpaired' }
    trim_failed              { path 'trimmed_reads/failed' } 
    trim_html                { path 'trimmed_reads/QC/fastp' }
    trim_json                { path 'trimmed_reads/QC/fastp' }
    multiqc_trim_report      { path 'summary_reports' }
    multiqc_trim_data        { path 'trimmed_reads/QC/multiqc' }
    // Alignment and QC
    bam_minfilt              { path 'alignments/minfilt' }
    bam_filt                 { path 'alignments/filtered' }
    qualimap_results          { path 'alignments/filtered/QC/qualimap' }
    multiqc_qualimap_report  { path 'summary_reports' }
    multiqc_qualimap_data    { path 'alignments/filtered/QC/multiqc' }

}