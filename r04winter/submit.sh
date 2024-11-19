#!/bin/bash
#SBATCH --time=168:00:00
#SBATCH --ntasks=4
#SBATCH --mem-per-cpu=2000
#SBATCH --mail-type=ALL
#SBATCH --mail-user=tylarmurray@usf.edu
#SBATCH --array=1-2
#
# slurm submit script

module load apps/singularity/3.5.0

singularity run ../run.sif winter $SLURM_ARRAY_TASK_ID