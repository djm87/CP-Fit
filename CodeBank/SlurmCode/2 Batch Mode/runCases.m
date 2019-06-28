function [] = runCases(lc)
%runCases Calls the windows or linux executable in parallel
    if and(ispc,lc.options{3,2})
       [~,~]=system('CallParallel.bat')
    elseif isunix & lc.options{3,2} & lc.options{6,2}
       [~,~]=system('sbatch --wait CallParallel.sh') 
    elseif and(isunix,lc.options{3,2})
       [~,~]=system('bash CallParallel.sh') 
    end
end

