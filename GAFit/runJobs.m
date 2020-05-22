function runJobs(isunix,slurmFlag,systemParams,nPop,totalCases,runFoldName,exeName)

if (isunix && (slurmFlag == 1))
    WriteCallParallelSlurm(systemParams,nPop);
    [~,~] = system('sbatch --wait CallParallelSlurm.sh');
    
elseif (slurmFlag == 2)
    jobList = cell(nPop*totalCases,1);
    jobCommand = cell(nPop*totalCases,1);
    jobCount = 1;
    for i = 1:nPop
        for j = 1:totalCases
            jobList{jobCount} = [runFoldName,'/',num2str(i),'/',num2str(j)];
            jobCommand{jobCount} = ['cd /d ',jobList{jobCount},' && ',exeName];
            jobCount = jobCount + 1;
        end
    end
    
    parfor i = 1:length(jobList)
        s = dir([jobList{i},'/STR_STR.OUT']);
        if (isempty(s) || s.bytes < 20000)
            [~,~] = system(jobCommand{i});
        end
    end
elseif (slurmFlag == 3)
    jobList = cell(nPop*totalCases,1);
    jobCommand = cell(nPop*totalCases,1);
    jobCount = 1;
    for i = 1:nPop
        for j = 1:totalCases
            jobList{jobCount} = [runFoldName,'/',num2str(i),'/',num2str(j)];
            jobCommand{jobCount} = ['cd ',jobList{i},' && ./',exeName];
            jobCount = jobCount + 1;
        end
    end
    
    parfor i = 1:length(jobList)
        [~,~] = system(jobCommand{i});
    end
end

end