#!/bin/bash

#SBATCH --job-name=WE43_Maud_Refinements
#SBATCH --partition=thrust2
#SBATCH --nodes=1
#SBATCH --tasks-per-node=7
#SBATCH --time=01:00:00             # Time limit hrs:min:sec
#SBATCH --output=WE43_Maud_Refinements_%A-%a.out    # Standard output and error log
#SBATCH -e array_%A-%a.err
#SBATCH --array=0-0                 # Array range

declare -a runStart=(1)
declare -a runEnd=(7)
seq ${runStart[${SLURM_ARRAY_TASK_ID}]} ${runEnd[${SLURM_ARRAY_TASK_ID}]} | xargs -n 1 -P 7 bash Maud_batch.sh
