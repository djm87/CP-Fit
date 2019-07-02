function WriteGAFitBatch(info,cases,GAoptions)

fid=fopen('CallParallel.sh','w');
fprintf(fid,'#! /bin/sh -f\n\n');
fprintf(fid,'#SBATCH --job-name=%s\n',runName);
fprintf(fid,'#SBATCH --partition=%s\n',info.sysInfo.slurmPartition);
fprintf(fid,'#SBATCH --nodes=1\n');
fprintf(fid,'#SBATCH --tasks-per-node=%1d\n',ncores);
% fprintf(fid,'#SBATCH --mem-per-cpu=%dgb           # Memory per processor\n',maxRam);
fprintf(fid,'#SBATCH --time=%s             # Time limit hrs:min:sec\n',maxTime);
fprintf(fid,'#SBATCH --output=%s_%%A-%%a.out    # Standard output and error log\n',runName);
fprintf(fid,'#SBATCH -e array_%%A-%%a.err\n');
fprintf(fid,'#SBATCH --array=0-%d                 # Array range\n\n',narrays-1);

% fprintf(fid,'find . -name \\*CPU* -type f -delete\n\n');
fprintf(fid,strcat('declare -a runStart=(',repmat('%d ',1,narrays),')\n'),startRun);
fprintf(fid,strcat('declare -a runEnd=(',repmat('%d ',1,narrays),')\n'),endRun);
fprintf(fid,'seq ${runStart[${SLURM_ARRAY_TASK_ID}]} ${runEnd[${SLURM_ARRAY_TASK_ID}]} | xargs -n 2 -P %d bash GARunModel_batch.sh\n',ncores);
fclose(fid);

end