#!/bin/bash
#SBATCH --time=168:00:00
#SBATCH --ntasks=32
#SBATCH --mem-per-cpu=2000
#
# slurm submit script

module load apps/singularity/3.5.0

singularity run run.sif 
