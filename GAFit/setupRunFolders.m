function setupRunFolders(cases,populationSize,runSource)

runFold = 'RunningFolder';
if (exist(runFold,'dir') == 7)
    disp('RunningFolder already exist.');
else
    mkdir(runFold);
end

caseIDs = cases{1:end,'CaseIdentifier'};

sourceDir = dir(runSource);
sourceDir([sourceDir.isdir]) = [];

% Make one folder for each population for each case specified by
% CaseIdentifier column in the WhatToFit file
for i = 1:length(caseIDs)
    curCase = [runSource,'/',caseIDs{i}];
    caseDir = dir(curCase);
    caseDir([caseDir.isdir]) = [];
    for j = 1:populationSize
        curDir = [runFold,'/',caseIDs{i},'/',num2str(j)];
        if (exist(curDir,'dir') ~= 7)
            mkdir(curDir);
        end
        % Copy the contents of the source files folder into each run folder
        % These are the files common to all runs
        for k = 1:numel(sourceDir)
            sourceFile = [sourceDir(k).folder,'/',sourceDir(k).name];
            destFile = [curDir,'/',sourceDir(k).name];
            copyfile(sourceFile,destFile);
        end
        
        % Copy the contents of the run files specific to each case from the
        % folders specified by CaseIdentifier column inside the source file
        % folder
        for k = 1:numel(caseDir)
            sourceFile = [caseDir(k).folder,'/',caseDir(k).name];
            destFile = [curDir,'/',caseDir(k).name];
            copyfile(sourceFile,destFile);
        end
    end
end

end