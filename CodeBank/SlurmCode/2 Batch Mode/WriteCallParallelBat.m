function [] = WriteCallParallelBat(lc)
%writeCallParallelBat To get around the terrible parralel management in
%windows, and poor parallel support in octave, a script is used to write each
%case of concurrent runs and in the case of linux xargs is used in parallel. 
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
  elseif isunix
    ncores=num2str(lc.options{1,2});
    fid=fopen('CallParallel.sh','w');
    fprintf(fid,'#!/bin/sh\n\n');
    fprintf(fid,'seq 1 %s | xargs -n 1 -P %s bash Maud_batch.sh\n\n',ncores,ncores);
    fclose(fid);      
  end
    

end