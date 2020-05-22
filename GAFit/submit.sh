#!/bin/bash
#
#SBATCH --job-name=runGA
#SBATCH --nodes=1
#SBATCH --tasks-per-node=4
##SBATCH --cpus-per-task=1
#SBATCH --partition=thrust2
#SBATCH -s
#SBATCH --mem=300G     # Memory
#SBATCH -o ./headlessOut/%A_%a.out
#SBATCH -e ./headlessOut/%A_%a.err
##SBATCH --array=1-16%4
#SBATCH -t 60-0:0:0
#SBATCH --time-min=1:00:00
#SBATCH --mail-user=zf1005@wildcats.unh.edu
#SBATCH --mail-type=START,END

#echo "My SLURM_ARRAY_TASK_ID: " $SLURM_ARRAY_TASK_ID

#Load the required modules for Abaqus with fortran
module load MATLAB

./runHeadless.sh MatlabRun
