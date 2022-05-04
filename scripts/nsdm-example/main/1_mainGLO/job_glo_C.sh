#!/bin/bash
#SBATCH --job-name=glo_C
#SBATCH --output=%x_%A-%a.out
#SBATCH --error=%x_%A-%a.err

Rscript glo_C.R $SLURM_ARRAY_TASK_ID