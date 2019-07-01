function errTotal = modelWrapErrFunction(par,cases,fitParam,info)

%% Global Variables
global runData; 
global runGeneration;
runGeneration = runGeneration + 1;
runData.params{runGeneration} = par;

%% Shared variable localization
fitRange = [cases.Start,cases.End];
Table = fitParam;
ishift = info.fitStrat.ishift;
isystem = info.sysInfo.sysType;
nPop = size(par,1);

%% Write sx files
caseIDs = cases{1:end,'CaseIdentifier'};
paramFname = info.modelInfo.paramFileName;

% write the sx files for each specific case then they will be copied into
% each population run folder
paramScaling = Table.scaling;
parfor i = 1:nPop
    for j = 1:numel(caseIDs)
        % writes the parameter file nPop times for each case
        WriteSxFile(['RunningFolder/',caseIDs{j},'/',num2str(i),'/',paramFname],[par(i,:)',paramScaling],caseIDs{j});
    end
end

%% Write slurm batch scripts and execute the scripts
% The batch script should run each executable in the RunningFolder

%% 
curSimData = cell(nPop*length(caseIDs),1);
activitiesPH1 = cell(nPop*length(caseIDs),1);
activitiesPH2 = cell(nPop*length(caseIDs),1);
activitiesPH3 = cell(nPop*length(caseIDs),1);
errors = zeros(nPop*length(caseIDs),1);
runData.ssCurves{runGeneration} = cell(nPop,length(caseIDs));
runData.err{runGeneration} = zeros(nPop,length(caseIDs));

outputFileName = cases.SimOut;
colX = cases.ColumnX;
colY = cases.ColumnY;
dataFiles = cases.FilePath;
fitRange = [cases.Start,cases.End];
parfor i = 1:nPop
    for j = 1:numel(caseIDs)
        % save the data
        vpscData = importdata(['RunningFolder/',caseIDs{j},'/',num2str(i),'/',outputFileName{j}]); %#ok<PFBNS>
        vpscData = vpscData.data(colX(j),colY(j)); %#ok<PFBNS>
        curSimData{i} = vpscData;
        
        % Calculate error ------------------------------------------------
        simModx = abs(vpscData(:,1));
        simMody = abs(vpscData(:,2));
        
        % Determine the case
        expData = importdata(dataFiles{j}); %#ok<PFBNS>
        expX = expData(:,1);
        expY = expData(:,2);
        
        % Fit simulated model data to a function
        [simModelx, simModely] = prepareCurveData(simModx,simMody);
        [simFitFcn, ~] = fit(simModelx,simModely,'smoothingspline');
        
        errors(i) = calcError(ishift,expX,expY,simFitFcn,fitRange(j)); %#ok<PFBNS>
        
        vpscData = importdata(['RunningFolder/',caseIDs{j},'/',num2str(i),'/ACT_PH1.OUT']);
        vpscData = vpscData.data;
        activitiesPH1{i} = vpscData;
        vpscData = importdata(['RunningFolder/',caseIDs{j},'/',num2str(i),'/ACT_PH2.OUT']);
        vpscData = vpscData.data;
        activitiesPH2{i} = vpscData;
        vpscData = importdata(['RunningFolder/',caseIDs{j},'/',num2str(i),'/ACT_PH3.OUT']);
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