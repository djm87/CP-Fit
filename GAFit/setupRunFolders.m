function setupRunFolders(cases,populationSize,runSource,runFoldName)

if (exist(runFoldName,'dir') ~= 7)
    mkdir(runFoldName);
end

caseIDs = cases{1:end,'CaseIdentifier'};

for j = 1:length(caseIDs)
    copySource = [runSource,'/',caseIDs{j},'/'];
    parfor i = 1:populationSize
        copyDestination = [runFoldName,'/',num2str(i),'/',num2str(j),'/'];
        if (exist(copyDestination,'dir') ~= 7)
            mkdir(copyDestination);
            copyfile(copySource,copyDestination);
        end
    end
end


end