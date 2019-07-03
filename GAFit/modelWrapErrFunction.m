function errTotal = modelWrapErrFunction(par,systemParams)

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
fitRange = [cases.Start,cases.End];
Table = fitParam;
ishift = info.fitStrat.ishift;
par
nPop = size(par,1)
totalCases = length(caseIDs);

%% Write sx files
% write the sx files for each specific case then they will be copied into
% each population run folder
paramScaling = Table.scaling;
input = Table.Parameters;
for i = 1:nPop
    i
    input(Table.fitFlag == true) = par(i,:)';
    for j = 1:numel(caseIDs)
        j
        % writes the parameter file nPop times for each case
        WriteSxFile([runFoldName,'/',num2str(i),'/',num2str(j),'/',paramFname],[input,paramScaling],caseIDs{j});
    end
end

%% Write slurm batch scripts and execute the scripts
% The batch script should run each executable in the RunningFolder
if (isunix && (info.sysInfo.slurmFlag == 1))
    WriteCallParallelSlurm(systemParams);
    [~,~] = system('sbatch --wait CallParallelSlurm.sh');
end

%% 
curSimData = cell(nPop*totalCases,1);
activitiesPH1 = cell(nPop*totalCases,1);
activitiesPH2 = cell(nPop*totalCases,1);
activitiesPH3 = cell(nPop*totalCases,1);
errors = zeros(nPop*totalCases,1);
runData.ssCurves{runGeneration} = cell(nPop,totalCases);
runData.err{runGeneration} = zeros(nPop,totalCases);

outputFileName = cases.SimOut;
colX = cases.ColumnX;
colY = cases.ColumnY;
dataFiles = cases.FilePath;
for i = 1:nPop
    for j = 1:numel(caseIDs)
        % save the data
        vpscData = importdata([runFoldName,'/',num2str(i),'/',num2str(j),'/',outputFileName{j}]);
        vpscData = [vpscData.data(:,colX(j)),vpscData.data(:,colY(j))];
        curSimData{i} = vpscData;
        
        % Calculate error ------------------------------------------------
        simModx = abs(vpscData(:,1));
        simMody = abs(vpscData(:,2));
        
        % Determine the case
        expData = importdata(dataFiles{j});
        expX = expData(:,1);
        expY = expData(:,2);
        
        % Fit simulated model data to a function
        [simModelx, simModely] = prepareCurveData(simModx,simMody);
        [simFitFcn, ~] = fit(simModelx,simModely,'smoothingspline');
        
        errors(i) = calcError(ishift,expX,expY,simFitFcn,fitRange(j,:),simModelx);
        
        vpscData = importdata([runFoldName,'/',num2str(i),'/',num2str(j),'/ACT_PH1.OUT']);
        vpscData = vpscData.data;
        activitiesPH1{i} = vpscData;
        vpscData = importdata([runFoldName,'/',num2str(i),'/',num2str(j),'/ACT_PH2.OUT']);
        vpscData = vpscData.data;
        activitiesPH2{i} = vpscData;
        vpscData = importdata([runFoldName,'/',num2str(i),'/',num2str(j),'/ACT_PH3.OUT']);
        vpscData = vpscData.data;
        activitiesPH3{i} = vpscData;
    end
end
runData.mainSimData{runGeneration} = reshape(curSimData,[nPop,totalCases]);
runData.activitiesPH1{runGeneration} = reshape(activitiesPH1,[nPop,totalCases]);
runData.activitiesPH2{runGeneration} = reshape(activitiesPH2,[nPop,totalCases]);
runData.activitiesPH3{runGeneration} = reshape(activitiesPH3,[nPop,totalCases]);
runData.err{runGeneration} = reshape(errors,[nPop,totalCases]);
errTotal = reshape(errors,[nPop,totalCases]);
errTotal = sum(errTotal,2);

end