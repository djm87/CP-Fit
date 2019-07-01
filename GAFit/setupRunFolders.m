function setupRunFolders(cases,populationSize,runSource)

runFold = 'RunningFolder';
if (exist(runFold,'dir') == 7)
    disp('RunningFolder already exist.');
else
    mkdir(runFold);
end

caseIDs = cases{1:end,'CaseIdentifier'};

% Make one folder for each population for each case specified by
% CaseIdentifier column in the WhatToFit file
for i = 1:length(caseIDs)
    curCase = [runSource,'/',caseIDs{i}];
    for j = 1:populationSize
        curDir = [runFold,'/',caseIDs{i},'/',populationSize];
        if (exist(curDir,'dir') ~= 7)
            mkdir(curDir);
        end
        % Copy the contents of the source files folder into each run folder
        % These are the files common to all runs
        copyfile(runSource,curDir);
        
        % Copy the contents of the run files specific to each case from the
        % folders specified by CaseIdentifier column inside the source file
        % folder
        copyfile(curCase,curDir);
    end
end

end