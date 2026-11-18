# kelpflow
nextflow pipeline for analysis of kelp whole genome data 
written by Meaghan I. Clark

This pipeline was designed to run quality control and basic filtering on low/mixed coverage whole genome sequencing data for kelp projects.

The major steps are as following: 

### 1. DATA PROCESSING AND QUALITY ASSESSMENT
        a. Runs FASTQC on raw reads and summarizes results into a multiQC report
        b. trims reads using fastp and summarizes results into a multiQC report
        c. aligns reads to a reference genome using bwa-mem2 and proceses alignments using samtools
        d. runs qualimap to assess alignment qualit and summarizes results into a multiQC report
### 2. FILTERING
        a. filters individuals based on a minimum average coverage cut-off
        b. runs ngsParalog to identify regions of problematic mapping 
        c. generates a depth histogram for visual assessment of coverage for secondary cut-off
        d. runs ngsLD and generates a list of putative variants sites that passed both mapping and LD filtering. 
### 3. ANALYSES
        a. generates a pca using pcangsd for initial assessment of population structure and/or batch effects

To use on elkhorn, use this nextflow orchestrator script: 

sbatch run_nextflow.sh

## SET UP
This pipeline uses Elkhorn's module system, conda environments, local installs, and singularity containers. Maybe someday in the future there will be a dedicated singularity container, but that day is not today.

### Required programs

These are the programs required by the pipeline. Nextflow will load everything for you, including creating conda environments. I have included the module information below for future troubleshooting purposes. 
 
- miniconda3 (elkhorn module, `module load miniconda3/3.13`)
- nextflow (elkhorn module, `module load nextflow/25.10.4`)
- R (elkhorn module, `module load r/4.4.1`)
- fastQC (elkhorn module, `module load fastqc/0.12.1`)
- multiqc (elkhorn module, `module load multiqc/1.27`)
- fastp (elkhorn module, `module load fastp/0.23.2`)
- bwa-mem2 (elkhorn module, `module load bwa-mem2/2.2.1`)
- samtools (elkhorn module, `module load samtools/1.21`)
- bamUtil (conda environment)
- qualimap (conda environment)
- angsd (elkhorn module, `module load angsd/0.940`)
- singularity (elkhorn module, `module load singularity-ce/singularity-ce.4.1.4`)
- ngsParalog (singularity container, see below*)
- bedtools (elkhorn module, `bedtools/2.26.0`)
- ngsLD (elkhorn module, `module load ngsld/1.2.0`)
- prune_graph (user install, see below*)
- pcangsd (conda environment. Also exists as an elkhorn module)

/* ngsParalog and prune_graph require some additional steps from the user before the pipeline is able to run. 

### ngsParalog installed as a singularity container
Singularity containers need to be set up before running the pipeline, using the code below. Run the following commands from the kelpflow directory. 

    #### Singularity container for ngsparalog
        ```
        module load singularity-ce/singularity-ce.4.1.4

        singularity build --fakeroot containers/ngsparalog.sif containers/ngsparalog.def
        singularilty build --sandbox containers/ngsparalog_sandbox containers/ngsparalog.sif
        ```

Then, change the `TMPDIR` and `SINGULARITY_TMPDIR` parameters in nextflow.config. 

***TMPDIR and SINGULARITY_TMPDIR in ngsparalog_train_hmm.nf deeds to be modified to the user's preference. Meaghan should make this a parameter in nextflow.config! 

### prune_graph needs to be installed locally: 

    #### Local installation of prunegraph
    # https://github.com/fgvieira/prune_graph

    ```
    module load git 
    git clone https://github.com/fgvieira/prune_graph.git
    ```

    After cloning the git repository, update the file path to your local install in nextflow.config

### HOW TO: TESTING
    nextflow run preview.nf # use preview.nf to view channels when testing 

### HOW TO: MAKING FLOWCHART
    module load nextflow/25.10.4

    nextflow run main.nf -preview -with-dag flowchart.mmd # change extension for different file type, i.e., ".png" 
    nextflow run main.nf -preview -with-dag flowchart.png 

    # edit using nf-metro-playground

### HOW TO: Creating new singularity containers

    # start interactive job: salloc --time=168:00:00 --mem=10Gb --ntasks=1 --cpus-per-task=1 --partition=windfall

    module load singularity-ce/singularity-ce.4.1.4

    singularity build --fakeroot containers/ngsparalog.sif containers/ngsparalog.def
    singularilty build --sandbox containers/ngsparalog_sandbox containers/ngsparalog.sif
    
    #test
    singularity exec ngsparalog.sif ngsParalog
    singularity exec ngsparalog.sif bash -c "source activate dupHMM && Rscript -e 'library(truncnorm); library(docopt); library(expm); library(codetools)'"

# NOTES AND FUTURE FEATURES
* add k-mer step to get genome profile as part of QC step? 
* add reference QC/evaluation step 

