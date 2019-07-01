#!/bin/bash
#SBATCH --job-name=WE43_Maud_refinements
#SBATCH --partition=thrust2
#SBATCH --nodes=1
#SBATCH --tasks-per-node=32
#SBATCH --mem-per-cpu=8gb           # Memory per processor
#SBATCH --time=00:30:00             # Time limit hrs:min:sec
#SBATCH --output=array_%A-%a.out    # Standard output and error log
#SBATCH -e array_%A-%a.err
#SBATCH --array=0-3                 # Array range

find . -name \*CPU* -type f -delete

declare -a runStart=(1 33 65 97)
declare -a runEnd=(32 64 96 106)
seq ${runStart[${SLURM_ARRAY_TASK_ID}]} ${runEnd[${SLURM_ARRAY_TASK_ID}]} | xargs -n 1 -P 32 bash Maud_batch.sh
#bash Maud_batch.sh ${SLURM_ARRAY_TASK_ID}
