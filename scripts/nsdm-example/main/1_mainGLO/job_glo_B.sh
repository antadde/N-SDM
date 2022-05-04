#!/bin/bash
#SBATCH --job-name=glo_B
#SBATCH --output=%x_%A-%a.out
#SBATCH --error=%x_%A-%a.err

export OMP_NUM_THREADS=1

Rscript glo_B.R $SLURM_ARRAY_TASK_ID