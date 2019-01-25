#!/bin/bash
#SBATCH --job-name s50 
#SBATCH --partition short 
#SBATCH --ntasks 1 
#SBATCH --time 0-12:00        # Wall time limit in D-HH:MM
#SBATCH --mem-per-cpu=3000     # Memory limit for each tasks (in MB)
#SBATCH -o out/Rt_%j.out    # File to which STDOUT will be written
#SBATCH -e out/Rt_%j.err    # File to which STDERR will be written

module load R/3.5.1-foss-2018b
Rscript simulation-50.R $SLURM_ARRAY_TASK_ID > ../experiments/log/g50-$SLURM_ARRAY_TASK_ID

# sbatch --array=1-400 simulation-50.sbatch
