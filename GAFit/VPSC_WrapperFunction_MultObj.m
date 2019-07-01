function errTotal = VPSC_WrapperFunction_MultObj(par,systemParams)
    
    %% Run write Sx File and immediately a.out in parallel
    parfor i = 1:size(targetFolders,1)

        % Move into the EPSC folder to run a.out or the windows exe
        cd(targetFolders{i});
        if (isystem == 0)
            % run the compression files in batch for windows
            [~,~] = system('epsc4_6_4_vol_avg_spin.exe');
        elseif (isystem == 1)
            % run the compression files in batch for linux
            [~,~] = system([curFold,'/a.out']);
        end
        
        
    end
end