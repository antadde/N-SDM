#!/bin/bash
#SBATCH --job-name=pre_A
#SBATCH --output=./0_mainPRE/%x_%A-%a.out
#SBATCH --error=./0_mainPRE/%x_%A-%a.err

export OMP_NUM_THREADS=1

Rscript ./0_mainPRE/pre_A.R