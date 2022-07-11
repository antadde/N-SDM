#!/bin/bash
#SBATCH --job-name=fut_B
#SBATCH --output=%x_%A-%a.out
#SBATCH --error=%x_%A-%a.err

Rscript fut_B.R $SLURM_ARRAY_TASK_ID