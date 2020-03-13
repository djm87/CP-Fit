function setupRunFolders(cases,populationSize,runSource,runFoldName)

if (exist(runFoldName,'dir') == 7)
    disp('RunningFolder already exist.');
else
    mkdir(runFoldName);
    caseIDs = cases{1:end,'CaseIdentifier'};
    
    copySource = cell(length(caseIDs,1),1);
    copyDestination = cell(populationSize);
    for i = 1:length(caseIDs)
        copySource{i} = [runSource,'/',caseIDs{i},'/'];
        for j = 1:populationSize
            copyDestination{j} = [runFoldName,'/',num2str(j),'/',num2str(i),'/'];
        end
    end
    
    for i = 1:length(caseIDs)
        curCopySource = copySource{i};
        parfor j = 1:populationSize
            mkdir(copyDestination{j});
            copyfile(curCopySource,copyDestination{j});
        end
    end
    
end

end