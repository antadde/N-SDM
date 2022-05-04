#!/bin/bash
#SBATCH --job-name=FUT_C
#SBATCH --output=%x_%A-%a.out
#SBATCH --error=%x_%A-%a.err

Rscript fut_C.R $SLURM_ARRAY_TASK_ID