nextflow.enable.dsl = 2


// Import modules
// set up
include { SPLIT_GENOME                } from './modules/local/split_genome.nf'
// data processing
include { FASTQC                      } from './modules/local/fastqc.nf'
include { MULTIQC as MULTIQC_RAW      } from './modules/local/multiqc.nf'
include { MULTIQC as MULTIQC_TRIM     } from './modules/local/multiqc.nf'
include { TRIM_READS as TRIM          } from './modules/local/trim_reads.nf'
include { ALIGN                       } from './modules/local/align_to_genome.nf'
include { QUALIMAP                    } from './modules/local/qualimap.nf'
include { MULTIQC as MULTIQC_QUALIMAP } from './modules/local/multiqc.nf'
include { PARSE_QUALIMAP              } from './modules/local/parse_qualimap.nf'
// filtering
include { FILTER_BY_COVERAGE          } from './modules/local/filter_by_coverage.nf'
include { MAKE_NGSPARALOG_INPUT       } from './modules/local/ngsparalog_make_input.nf'
include { CALC_LR                     } from './modules/local/ngsparalog_calc_LR.nf'
include { TRAIN_HMM                   } from './modules/local/ngsparalog_train_hmm.nf'
include { RUN_HMM                     } from './modules/local/ngsparalog_run_hmm.nf'
include { DEPTH_HIST                  } from './modules/local/depth_hist.nf'
include { GENERATE_GENO_LIKELIHOODS   } from './modules/local/generate_geno_likelihoods.nf'
include { FILTER_PARALOGS             } from './modules/local/filter_paralogs.nf'
include { RUN_NGSLD                   } from './modules/local/run_ngsld.nf'
include { PRUNE_LD_GRAPH              } from './modules/local/prune_graph.nf'
// Analyses
include { GENERATE_PCA_LIKELIHOODS    } from './modules/local/generate_pca_likelihoods.nf'
include { RUN_PCANGSD                 } from './modules/local/run_pcangsd.nf'


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
    // build reference index channel
    ch_fai = ch_ref_genome
        .map { meta, fasta, index_only -> index_only.find { it.name.endsWith('.fai') } }


    // Start pipeline

    // Set up: split genome into 20 segments for parallel processing 
    SPLIT_GENOME(ch_fai)

    // CHANNEL: Segments of the genome after splitting
    ch_segments = SPLIT_GENOME.out.segments
        .flatten()
        .map { f -> tuple(f.baseName.replaceAll('segment_', ''), f) }
    
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

    // CHANNEL: Gather minimally filtered alignment files (for ngsParalog)
    ch_minfilt_bams = ALIGN.out.bam_minfilt

    // Run Qualimap on filtered BAM files
    QUALIMAP(ch_filtered_bams.map { meta, bam, bai -> tuple(meta, bam) })

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

    // Filter samples by coverage
    FILTER_BY_COVERAGE(ch_qualimap_stats)

    // CHANNEL: individuals that passed coverage filter
    ch_passing_inds = FILTER_BY_COVERAGE.out.pass_cov_samples
        .splitText()
        .map{ it.trim() }

    // CHANNEL: minfilt bams for passing individuals
    ch_passing_minfilt_bams = ALIGN.out.bam_minfilt
        .map { meta, bam, bai -> tuple(meta.id, bam, bai) }
        .combine(ch_passing_inds.map { id -> tuple(id, true) }, by: 0)
        .map { id, bam, bai, flag -> tuple(bam, bai) }
    ch_all_passing_minfilt_bams = ch_passing_minfilt_bams.map { bam, bai -> bam }.collect()     // Collect every passing sample's BAM/BAI into single lists (staged once, reused per segment)
    ch_all_passing_minfilt_bais = ch_passing_minfilt_bams.map { bam, bai -> bai }.collect()
    ch_minfilt_bamlist = ch_passing_minfilt_bams
        .map { bam, bai -> bam.name }
        .collectFile(name: 'bamlist_minfilt.txt', newLine: true)
        .first()

    // CHANNEL: filtered bams for passing individuals
    ch_passing_filt_bams = ALIGN.out.bam_filt
        .map { meta, bam, bai -> tuple(meta.id, bam, bai) }
        .combine(ch_passing_inds.map { id -> tuple(id, true) }, by: 0)
        .map { id, bam, bai, flag -> tuple(bam, bai) }
    ch_all_passing_filt_bams = ch_passing_filt_bams.map { bam, bai -> bam }.collect()     // Collect every passing sample's BAM/BAI into single lists (staged once, reused per segment)
    ch_all_passing_filt_bais = ch_passing_filt_bams.map { bam, bai -> bai }.collect()
    ch_filt_bamlist = ch_passing_filt_bams
        .map { bam, bai -> bam.name }
        .collectFile(name: 'bamlist_filt.txt', newLine: true)
        .first()

    // START: ngsParalog filter
    
    // make input files for ngsParalog
    MAKE_NGSPARALOG_INPUT(ch_segments, ch_minfilt_bamlist, ch_all_passing_minfilt_bams, ch_all_passing_minfilt_bais)

    // CHANNEL: gather ngsparalog input with genome segments
    ch_segment_posfiles = ch_segments
        .join(MAKE_NGSPARALOG_INPUT.out.varsites)

    // CHANNEL: package reference and index into tuple
    ch_ref_fasta_fai = ch_ref_genome
        .map { meta, fasta, index_only ->
            def fai = index_only.find { it.name.endsWith('.fai') }
            tuple(fasta, fai)
    }
    // Calculate ngsParalog likelihood ratios
    CALC_LR(
        ch_segment_posfiles,
        ch_minfilt_bamlist,
        ch_all_passing_minfilt_bams,
        ch_all_passing_minfilt_bais,
        ch_ref_fasta_fai
    )

    // CHANNEL: Combine all segments' LR values into one genome-wide file
    ch_lr_all = CALC_LR.out.lr
        .map { segment_id, tsv -> tsv }
        .collectFile(name: 'Mpyr_LR_all.txt', storeDir: "${params.outdir}/ngsparalog/LR")

    // CHANNEL: Combine all segments' coverage values into one genome-wide file
    ch_coverage_all = MAKE_NGSPARALOG_INPUT.out.coverage
        .map { segment_id, cov -> cov }
        .collectFile(name: 'average_ind_coverage_all.txt', storeDir: "${params.outdir}/ngsparalog/input")

    // Train HMM to estimate parameters 
    TRAIN_HMM(ch_lr_all, ch_coverage_all)

    // CHANNEL: Per-segment LR + coverage, paired by segment_id
    ch_segment_lr_cov = CALC_LR.out.lr
        .join(MAKE_NGSPARALOG_INPUT.out.coverage)
        // -> tuple(segment_id, lr_file, cov_file)
    
    // Run HMM using estimated parameters 
    RUN_HMM(ch_segment_lr_cov, TRAIN_HMM.out.paramfile.first())

    // CHANNEL: Combine HMM regions (bad regions) from all segements
    ch_all_paralogous_regions = RUN_HMM.out.hmm_regions
        .map { segment_id, regions_file -> regions_file }
        .collect()
    
    // CHANNEL: total depth
    ch_coverage_total = MAKE_NGSPARALOG_INPUT.out.coverage_total
        .map { segment_id, cov -> cov }
        .collectFile(name: 'total_ind_coverage_all.txt')

    // Depth histogram to inform filters
    DEPTH_HIST(ch_coverage_total)

    // LD pruning

    // CHANNEL: collect putative variable sites (minimal quality filtering) into position file
    ch_putative_sites = MAKE_NGSPARALOG_INPUT.out.varsites
        .map { segment_id, f -> f }
        .collectFile(name: 'putative_varsites_all.pos')

    // Remove areas with problematic mapping as identified by ngsParalog from list of putative variants
    FILTER_PARALOGS(ch_putative_sites, ch_all_paralogous_regions)

    // Generate genotype likelihoods for all passing putative variants 
    GENERATE_GENO_LIKELIHOODS(
        ch_filt_bamlist, ch_all_passing_filt_bams, ch_all_passing_filt_bais,
        FILTER_PARALOGS.out.sites, DEPTH_HIST.out.cutoff
    )

    // Run ngsLD to calculate pairwise LD between all passing putative variants
    RUN_NGSLD(GENERATE_GENO_LIKELIHOODS.out.ld_input)       

    // Run prune graph to remove linked sites from the list of putative variants
    PRUNE_LD_GRAPH(RUN_NGSLD.out.ld)

    // Make PCA

    // CHANNEL: pruned sites as a plain path
    ch_pruned_sites = PRUNE_LD_GRAPH.out.pruned_sites.map { meta, f -> f }

    // generate genotype likelihoods for PASS sites
    GENERATE_PCA_LIKELIHOODS(
        ch_filt_bamlist, ch_all_passing_filt_bams, ch_all_passing_filt_bais,
        ch_pruned_sites
    )

    // run pcangsd 
    RUN_PCANGSD(GENERATE_PCA_LIKELIHOODS.out.beagle)

    // Major outputs: 
        // list of PASS individuals
        // list of PASS sites
        // report summarizing data quality, alignment quality, filtering

    // Additional steps to consider adding or making their own pipelines: 
        // ngsRelate to filter related individuals 
        // admixture analyses
        // pop gen summary stats
        // selection scans implemented in angsd
        // GEA


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
    // Filtering
        samples_passing_coverage = FILTER_BY_COVERAGE.out.pass_cov_samples
        samples_failing_coverage = FILTER_BY_COVERAGE.out.fail_cov_samples
        depth_histogram          = DEPTH_HIST.out.histogram
        depth_cutoff             = DEPTH_HIST.out.cutoff
        hmm_regions              = RUN_HMM.out.hmm_regions
    // LD pruning
        beagle_likelihoods       = GENERATE_GENO_LIKELIHOODS.out.beagle
        ngsld_output             = RUN_NGSLD.out.ld
        ld_pruning_report        = PRUNE_LD_GRAPH.out.report
        pruned_safe_sites = PRUNE_LD_GRAPH.out.pruned_sites
    // PCA
        pca_beagle  = GENERATE_PCA_LIKELIHOODS.out.beagle
        pca_results = RUN_PCANGSD.out.all
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
    qualimap_results         { path 'alignments/filtered/QC/qualimap' }
    multiqc_qualimap_report  { path 'summary_reports' }
    multiqc_qualimap_data    { path 'alignments/filtered/QC/multiqc' }
    // Filtering
    samples_passing_coverage { path 'filters' }
    samples_failing_coverage { path 'filters' }
    depth_histogram          { path 'filters' }
    depth_cutoff             { path 'filters' }
    hmm_regions              { path 'filters' }
    // LD pruning
    beagle_likelihoods       { path 'ld_pruning/beagle' }
    ngsld_output             { path 'ld_pruning/ngsLD' }
    ld_pruning_report        { path 'ld_pruning/report' }
    pruned_safe_sites        { path 'ld_pruning/final_sites' }
    // PCA 
    pca_beagle               { path 'pca/beagle' }
    pca_results              { path 'pca/results' }

}