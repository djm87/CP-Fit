#!/bin/bash
#SBATCH --nodes=2
#SBATCH --job-name=hellowrld2 
#SBATCH --tasks-per-node=32
#SBATCH -o hellowrld2.out
#SBATCH --time=60-0:0:0
#SBATCH --time-min=10:00:00
#SBATCH --partition=thrust2

matlab -nodesktop -nosplash -nodisplay -r 'run helloworld.m; quit'