#!/bin/bash

#SBATCH --job-name=SlurmTest
#SBATCH --partition=thrust2
#SBATCH --nodes=1
#SBATCH --tasks-per-node=7
#SBATCH --time=01:00:00             # Time limit hrs:min:sec
#SBATCH --output=SlurmTest.out    # Standard output and error log
#SBATCH -e array_Test.err
#SBATCH --array=0-0                 # Array range

declare -a runStart=(1)
declare -a runEnd=(7)
seq ${runStart[${SLURM_ARRAY_TASK_ID}]} ${runEnd[${SLURM_ARRAY_TASK_ID}]} | xargs -n 2 -P 0 bash GAFit_batch.sh
