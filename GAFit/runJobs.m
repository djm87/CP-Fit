function runJobs(isunix,slurmFlag,systemParams,nPop,totalCases,runFoldName,exeName)

if (isunix && (slurmFlag == 1))
    WriteCallParallelSlurm(systemParams,nPop);
%     disp('here after write parallel');
    [~,~] = system('sbatch --wait CallParallelSlurm.sh');
%     disp('here after finished running parallel slurm');
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
        [~,~] = system(jobCommand{i});
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