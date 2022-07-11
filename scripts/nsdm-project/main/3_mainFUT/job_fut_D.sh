#!/bin/bash
#SBATCH --job-name=fut_D
#SBATCH --output=%x_%A-%a.out
#SBATCH --error=%x_%A-%a.err

Rscript fut_D.R $SLURM_ARRAY_TASK_ID