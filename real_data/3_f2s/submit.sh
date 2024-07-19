#!/bin/bash

#SBATCH --array=1-88
#SBATCH --mem=31G
#SBATCH --time=72:00:00
#SBATCH --partition=ncpu

i=$((${SLURM_ARRAY_TASK_ID}-1))

Rscript ./calc.R $i
Rscript ./calc_fine.R $i

