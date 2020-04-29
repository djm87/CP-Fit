function WriteCallParallelSlurm(systemParams,nPop)
%writeCallParallelBat To get around the terrible parralel management in
%windows, and poor parallel support in octave, a script is used to write each
%case of concurrent runs and in the case of linux xargs is used in parallel
%mode. %Support has also been added for slurm work jobs.
cases = systemParams{1};
info = systemParams{4};
% GAoptions = systemParams{5};
runFoldName = systemParams{6};

% nPop = GAoptions.PopulationSize;
nCase = length(cases.CaseIdentifier);
ncores = info.sysInfo.parTasks;
narrays = ceil(nPop*nCase/ncores);
maxRam = info.sysInfo.slurmMaxRam;
maxTime = info.sysInfo.slurmMaxTime;
runName = info.sysInfo.slurmJobName;
exeFname = info.sysInfo.exeName;
nParJobs = info.sysInfo.slurmNParallelJobs;

fid=fopen('CallParallelSlurm.sh','w');
fprintf(fid,'#!/bin/bash\n');
fprintf(fid,'#SBATCH --job-name=%s\n',runName);
fprintf(fid,'#SBATCH --partition=%s\n',info.sysInfo.slurmPartition);
fprintf(fid,'#SBATCH --nodes=1\n');
fprintf(fid,'#SBATCH --tasks-per-node=%1d\n',ncores);
fprintf(fid,'#SBATCH --mem=%dG     # Memory \n',maxRam);
fprintf(fid,'#SBATCH --time=%s                # Time limit hrs:min:sec\n',maxTime);
fprintf(fid,'#SBATCH -o /dev/null  # suppress output log\n');
fprintf(fid,'##SBATCH -o ./outputs/%s_%%A-%%a.out  # Standard output log\n',runName);
fprintf(fid,'##SBATCH -e ./outputs/%s_%%A-%%a.err        # Standard error log\n',runName);
fprintf(fid,'#SBATCH --array=0-%d%%%d           # Array range\n\n',narrays-1,nParJobs);

% fprintf(fid,'find . -name \\*CPU* -type f -delete\n\n');
fprintf(fid,'nFolders=%d\n',nPop);
fprintf(fid,'nCases=%d\n\n',nCase);

fprintf(fid,'for ((i=0; i<nFolders; i++));do\n');
fprintf(fid,'  jname=$((${i} + 1))\n');
fprintf(fid,'  for ((j=0; j<nCases; j++));do\n');
fprintf(fid,'     linIndx=$((${i} * ${nCases} + ${j}))\n');
fprintf(fid,'     cname=$((${j} + 1))\n');
fprintf(fid,'     JobNames[${linIndx}]="%s/%s/$jname/$cname/"\n',pwd,runFoldName);
fprintf(fid,'  done\n');
fprintf(fid,'done\n\n');

fprintf(fid,'jobLength=${#JobNames[@]}\n\n');

fprintf(fid,'beginning=$((${SLURM_ARRAY_TASK_ID} * ${SLURM_NTASKS_PER_NODE}))\n');
fprintf(fid,'ending=$((${SLURM_ARRAY_TASK_ID} * ${SLURM_NTASKS_PER_NODE} + ${SLURM_NTASKS_PER_NODE}))\n\n');

fprintf(fid,'if [ ${ending} -gt ${jobLength} ];then\n');
fprintf(fid,'ending=${jobLength}\n');
fprintf(fid,'fi\n\n');

fprintf(fid,'for ((i=beginning; i<ending; i++));do\n');
fprintf(fid,'   JobsOnThisNode[$i]=${JobNames[$i]}\n');
fprintf(fid,'done\n\n');

fprintf(fid,'for i in ${JobsOnThisNode[@]};do echo $i;done | xargs -n1 -P ${SLURM_NTASKS_PER_NODE} bash runExe.sh\n');
fclose(fid);

end