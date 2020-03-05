function errors = modelWrapErrFunction(par,systemParams)

%% Global Variables
global runData; 
global runGeneration;
runGeneration = runGeneration + 1;
runData.params{runGeneration} = par;

%% Shared variable localization
cases = systemParams{1};
fitParam = systemParams{2};
info = systemParams{3};
runFoldName = systemParams{5};

caseIDs = cases{1:end,'CaseIdentifier'};
paramFname = info.modelInfo.paramFileName;
Table = fitParam;
nPop = size(par,1);
totalCases = length(caseIDs);

%% Write sx files
% write the sx files for each specific case then they will be copied into
% each population run folder
paramScaling = Table.scaling;
input = Table.Parameters;
for i = 1:nPop
    input(Table.fitFlag == true) = par(i,:)';
    for j = 1:numel(caseIDs)
        % writes the parameter file nPop times for each case
        writeDPSxFile([runFoldName,'/',num2str(i),'/',num2str(j),'/'],paramFname,[input,paramScaling]);
    end
end

%% Write slurm batch scripts and execute the scripts
% The batch script should run each executable in the RunningFolder
runJobs(isunix,info.sysInfo.slurmFlag,systemParams,nPop,totalCases,runFoldName,info.sysInfo.exeName);

%% Call errorEvalWrap that calls various functions to evaluate error for each objective
curSimData = cell(nPop*totalCases,1);
% activitiesPH1 = cell(nPop*totalCases,1);
% activitiesPH2 = cell(nPop*totalCases,1);
% activitiesPH3 = cell(nPop*totalCases,1);
errors = zeros(nPop,totalCases);
runData.ssCurves{runGeneration} = cell(nPop,totalCases);
runData.err{runGeneration} = zeros(nPop,totalCases);

[curSimData,errors] = errorEvalWrap(cases,nPop,caseIDs,runFoldName,curSimData,info.cyclicFits,errors,info.fitStrat.ishift);

runData.mainSimData{runGeneration} = reshape(curSimData,[nPop,totalCases]);
% runData.activitiesPH1{runGeneration} = reshape(activitiesPH1,[nPop,totalCases]);
% runData.activitiesPH2{runGeneration} = reshape(activitiesPH2,[nPop,totalCases]);
% runData.activitiesPH3{runGeneration} = reshape(activitiesPH3,[nPop,totalCases]);
runData.err{runGeneration} = errors

% Sum of objective errors. Obsolete for multi-objective
% errTotal = reshape(errors,[nPop,totalCases]);
% errTotal = sum(errTotal,2);

end