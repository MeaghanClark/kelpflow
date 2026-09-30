# kelpflow
nextflow pipeline for analysis of kelp whole genome data 
written by Meaghan I. Clark

To use on elkhorn, use this nextflow orchestrator script: 

sbatch run_nextflow.sh

### TESTING
nextflow run preview.nf # use preview .nf to view channels when testing 

### Creating new singularity containers

    # start interactive job: salloc --time=168:00:00 --mem=10Gb --ntasks=1 --cpus-per-task=1 --partition=windfall
    
module load singularity-ce/singularity-ce.4.1.4

singularity build ngsparalog.sif ngsparalog.def

### MAKING FLOWCHART
nextflow run main.nf -preview -with-dag flowchart.mmd # change extension for different file type, i.e., ".png" 

    # edit using nf-metro-playground