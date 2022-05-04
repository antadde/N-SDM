#!/bin/bash
#SBATCH --job-name=end_A
#SBATCH --output=%x_%A-%a.out
#SBATCH --error=%x_%A-%a.err

export OMP_NUM_THREADS=1

Rscript end_A.R $SLURM_ARRAY_TASK_ID