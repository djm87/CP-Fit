function [lc] = SetCPUs(lc)
%SetCPUs If the numcores on machine < user specified cores, set to
%max machine cores. If ncases < cores set, set machine cores to ncases

        if ispc 
            [status numprocs]=system('echo %number_of_processors%');
        elseif isunix
            [status,numprocs] = system('nproc');
        end

        if and(lc.options{1,2}>numprocs,~lc.options{6,2})
            lc.options{1,2}=str2double(numprocs);
            X=sprintf('More concurrent cases specified than CPUs... setting option ''%s'' to max of %d',lc.options{1,1},numprocs);
            disp(X);
        end

        if lc.ncases{lc.BS}<lc.options{1,2} 
            lc.options{1,2} =lc.ncases{lc.BS}; 
        end
end