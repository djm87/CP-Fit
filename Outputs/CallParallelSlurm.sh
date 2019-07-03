#!/bin/bash
#SBATCH --job-name=GAFit
#SBATCH --partition=general
#SBATCH --nodes=1
#SBATCH --tasks-per-node=32
#SBATCH --time=1-00:00:00                # Time limit hrs:min:sec
#SBATCH --output=GAFit_%A-%a.out  # Standard output log
#SBATCH -e GAFit_%A-%a.err        # Standard error log
#SBATCH --array=0-04           # Array range

nFolders=3
nCases=3

for ((i=0; i<nFolders; i++));do
  jname=$((${i} + 1))
  for ((j=0; j<nCases; j++));do
     linIndx=$((${i} * ${nCases} + ${j}))
     cname=$((${j} + 1))
     JobNames[${linIndx}]="cd /mnt/home/thrust2/zf1005/Matlab/GAFit/RunningFolder/$jname/$cname/ && ./a.out"
  done
done

jobLength=${#JobNames[@]}

beginning=$((${SLURM_ARRAY_TASK_ID} * ${SLURM_NTASKS_PER_NODE}))
ending=$((${SLURM_ARRAY_TASK_ID} * ${SLURM_NTASKS_PER_NODE} + ${SLURM_NTASKS_PER_NODE}))

if [ ${ending} -gt ${jobLength} ];then
ending=${jobLength}
fi

for ((i=beginning; i<ending; i++));do
   JobsOnThisNode[$i]=${JobNames[$i]}
done

for i in ${JobsOnThisNode[@]};do echo $i;done | xargs -n1 -P ${SLURM_NTASKS_PER_NODE} bash runExe.sh
