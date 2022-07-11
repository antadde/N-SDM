#!/bin/bash
#SBATCH --job-name=glo_C
#SBATCH --output=%x_%A-%a.out
#SBATCH --error=%x_%A-%a.err

export OMP_NUM_THREADS=1

Rscript glo_C.R $SLURM_ARRAY_TASK_ID