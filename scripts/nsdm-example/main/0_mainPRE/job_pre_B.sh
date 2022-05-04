#!/bin/bash
#SBATCH --job-name=pre_B
#SBATCH --output=%x_%A-%a.out
#SBATCH --error=%x_%A-%a.err

export OMP_NUM_THREADS=1

Rscript pre_B.R