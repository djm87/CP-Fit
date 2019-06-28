function [] = WriteCallParallel(lc)
%writeCallParallelBat To get around the terrible parralel management in
%windows, and poor parallel support in octave, a script is used to write each
%case of concurrent runs and in the case of linux xargs is used in parallel 
%mode. %Support has also been added for slurm work jobs.
  if ispc
    fid=fopen('CallParallel.bat','w');
    fprintf(fid,'@echo off\n\n');
    fprintf(fid,'(\n');
    for i=1:length(lc.batName) 
      if lc.options{4,2}
        fprintf(fid,'start "MAUD_Batch_Instance" /b  call %s\n',lc.batName{i});
      else
        fprintf(fid,'start "MAUD_Batch_Instance" call %s\n',lc.batName{i}); 
      end
    end
    fprintf(fid,') | set /P "="\n\n');
    fclose(fid);
    
  elseif and(isunix,lc.options{6,2})
    ncores=lc.options{1,2};
    narrays=ceil(lc.ncases{lc.BS}/lc.options{1,2});
    maxRam=lc.options{7,2};
    maxTime=lc.options{9,2};
    runName=lc.options{10,2}
    startRun=1:ncores:lc.ncases{lc.BS};
    endRun=ncores:ncores:lc.ncases{lc.BS};
    if length(endRun)<narrays
        endRun(narrays)=lc.ncases{lc.BS}
    end
    
    
    fid=fopen('CallParallel.sh','w');
    fprintf(fid,'#!/bin/bash\n\n');
    fprintf(fid,'#SBATCH --job-name=%s\n',runName);
    fprintf(fid,'#SBATCH --partition=%s\n',lc.options{11,2});
    fprintf(fid,'#SBATCH --nodes=1\n');
    fprintf(fid,'#SBATCH --tasks-per-node=%1d\n',ncores);
%     fprintf(fid,'#SBATCH --mem-per-cpu=%dgb           # Memory per processor\n',maxRam);
    fprintf(fid,'#SBATCH --time=%s             # Time limit hrs:min:sec\n',maxTime);
    fprintf(fid,'#SBATCH --output=%s_%%A-%%a.out    # Standard output and error log\n',runName);
    fprintf(fid,'#SBATCH -e array_%%A-%%a.err\n');
    fprintf(fid,'#SBATCH --array=0-%d                 # Array range\n\n',narrays-1);

%     fprintf(fid,'find . -name \\*CPU* -type f -delete\n\n');
    fprintf(fid,strcat('declare -a runStart=(',repmat('%d ',1,narrays),')\n'),startRun);
    fprintf(fid,strcat('declare -a runEnd=(',repmat('%d ',1,narrays),')\n'),endRun);
    fprintf(fid,'seq ${runStart[${SLURM_ARRAY_TASK_ID}]} ${runEnd[${SLURM_ARRAY_TASK_ID}]} | xargs -n 1 -P %d bash Maud_batch.sh\n',ncores);
    fclose(fid); 
  elseif isunix
    ncores=num2str(lc.options{1,2});
    fid=fopen('CallParallel.sh','w');
    fprintf(fid,'#!/bin/sh\n\n');
    fprintf(fid,'seq 1 %s | xargs -n 1 -P %s bash Maud_batch.sh\n\n',ncores,ncores);
    fclose(fid);     
  end
end