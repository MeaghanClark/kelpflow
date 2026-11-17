# kelpflow
nextflow pipeline for analysis of kelp whole genome data 
written by Meaghan I. Clark

This pipeline was designed to run quality control and basic filtering on low/mixed coverage whole genome sequencing data for kelp projects.

To use on elkhorn, use this nextflow orchestrator script: 

sbatch run_nextflow.sh

## SET UP
This pipeline uses Elkhorn's module system, conda environments, local installs, and singularity containers. Singularity containers need to be set up before running the pipeline, using the code below. 

    #### Singularity container for ngsparalog
        ```
        module load singularity-ce/singularity-ce.4.1.4

        singularity build --fakeroot containers/ngsparalog.sif containers/ngsparalog.def
        singularilty build --sandbox containers/ngsparalog_sandbox containers/ngsparalog.sif
        ```
    
    #### Local installation of prunegraph
    # https://github.com/fgvieira/prune_graph

    ```
    module load git 
    git clone https://github.com/fgvieira/prune_graph.git
    ```

    After cloning the git repository, update the file path in nextflow.config

### HOW TO: TESTING
nextflow run preview.nf # use preview .nf to view channels when testing 

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

# NOTES
add k-mer step to get genome profile as part of QC step? 
add reference evaluation step 

