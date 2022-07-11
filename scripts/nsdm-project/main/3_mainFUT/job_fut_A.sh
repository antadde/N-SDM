#!/bin/bash
#SBATCH --job-name=fut_A
#SBATCH --output=%x_%A-%a.out
#SBATCH --error=%x_%A-%a.err

Rscript fut_A.R $SLURM_ARRAY_TASK_ID