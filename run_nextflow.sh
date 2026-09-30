#!/bin/bash
#SBATCH --job-name=nf_orchestrator
#SBATCH --time=24:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --output ./run_logs/nf_orchestrator-%j.log
#SBATCH --partition=lab-mpinsky
#SBATCH --qos=pi-mpinsky
#SBATCH --account=pi-mpinsky

module load nextflow/25.10.4

nextflow run main.nf -profile slurm --input_dir /home/miclark/kelpflow/test_data/raw --report_name MpyrAdapt_QC -resume