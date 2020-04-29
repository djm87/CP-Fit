function errors = modelWrapErrFunction(par,systemParams)

%% Global Variables
global runData; 
global runGeneration;
runGeneration = runGeneration + 1;
runData.params{runGeneration} = par;

%% Shared variable localization
cases = systemParams{1};
fitParam = systemParams{2};
fitRecipe = systemParams{3};
info = systemParams{4};
runFoldName = systemParams{6};
PFFits = systemParams{7};
expData = systemParams{8};

caseIDs = cases{1:end,'CaseIdentifier'};
paramFname = info.modelInfo.paramFileName;
nPop = size(par,1);
totalCases = length(caseIDs);
totalPF = sum(cases{1:end,'IsPF'});
totalObjectives = totalCases + totalPF;

%% Write sx files
% write the sx files for each specific case then they will be copied into
% each population run folder
paramScaling = fitParam.scaling;
parfor i = 1:nPop
    params = fitParam.Parameters;
    % remake the total list of parameters accounting for non-fitting terms
    params(fitRecipe) = par(i,:);
    
    for j = 1:numel(caseIDs)
        % writes the parameter file(s) of all cases for each population
        % writeDPSxFile([runFoldName,'/',num2str(i),'/',num2str(j),'/'],paramFname,[par(i,:)',paramScaling]);
        writeTiSxFile([runFoldName,'/',num2str(i),'/',num2str(j),'/'],paramFname,params.*paramScaling,cases{j,'Dataset'});
    end
end

%% Write slurm batch scripts and execute the scripts
% The batch script should run each executable in the RunningFolder
runJobs(isunix,info.sysInfo.slurmFlag,systemParams,nPop,totalCases,runFoldName,info.sysInfo.exeName);

%% Call errorEvalWrap that calls various functions to evaluate error for each objective
curSimData = cell(nPop,totalCases);
errors = zeros(nPop,totalCases);
shiftind = zeros(nPop,totalCases);
runData.err{runGeneration} = zeros(nPop,totalObjectives);

[curSimData,errors,shiftind] = errorEvalWrap(cases,nPop,caseIDs,runFoldName,curSimData,errors,shiftind,PFFits,expData,info.cyclicFits,info.fitStrat.ishift);

runData.err{runGeneration} = errors;
[~,lowestErr] = min(sum(runData.err{runGeneration,1},2));
runData.lowestErrSimData{runGeneration} = curSimData(lowestErr,:);
runData.shiftinds{runGeneration} = shiftind(lowestErr,:);

end