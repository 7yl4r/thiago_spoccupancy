#!/bin/bash
#SBATCH --time=168:00:00
#SBATCH --ntasks=4
#SBATCH --mem-per-cpu=8000
#SBATCH --array=1-2
#SBATCH --mail-type=END,FAIL,TIME_LIMIT
#SBATCH --mail-user=user@usf.edu
#
# slurm submit script

module load apps/singularity/3.5.0

singularity run ../run.sif winter $SLURM_ARRAY_TASK_ID
