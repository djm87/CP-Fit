function errTotal = VPSC_WrapperFunction_MultObj(par,systemParams)
    % This function evaluates the error between model and data
    % The wrapper interfaces with stress strain output from EPSC
    % The objective is run par which is a vecorized input to the function
    % on as many cores as possible simultaneously.

    %% Global Variables
    global runData; 
    global runGeneration;
%     global runType;
    runGeneration = runGeneration + 1;
    runData.params{runGeneration} = par;
    
    %% Shared variables breakdown
    % models = systemParams{1};
    rootPath = systemParams{2};
    modelPaths = systemParams{3};
    cases = systemParams{4};
    % nPop = systemParams{5};
    modelCaseDist = systemParams{6};
    expModNames = systemParams{7};
    fitRange = systemParams{8};
    iopt = systemParams{9};
    toRefine = systemParams{10};
    inputs = systemParams{11};
    % additional option items, see OptimizationStart.m for explainations
    ishift = iopt(1);
    ifitpart = iopt(2);
    isystem = iopt(3);
    
    %% Prepare all the folders, sx files and run all the a.out in parallel
    % Given par and known cases and models
    nPop = size(par,1);
    systemParams{5} = nPop;
    % function sets up the sub folders for each synchronous run
    setupFolders(systemParams);
    
    totalCases = sum(modelCaseDist);
    % Create a list of all the directories to be looped through, run
    % writesx and run a.out
    targetFolders = cell(nPop*totalCases,1);
    linLoading = cell(nPop*totalCases,1);
    linCases = cell(nPop*totalCases,1);
    foldCount = 0;
    for i = 1:size(modelPaths,1)
        for j = 1:size(cases{i},1)
            for k = 1:nPop
                foldCount = foldCount + 1;
                targetFolders{foldCount} = [modelPaths{i},'/',cases{i}(j,:),'_',int2str(k)];
                linCases{foldCount}=i;
                linLoading{foldCount}=cases{i}(j,:);
                
            end
        end
    end
    % Save the list of folders for deleting
    save([rootPath,'/del_fold.mat'],'targetFolders');
    
    curSSCurves = cell(nPop*totalCases,1);
    activitiesPH1 = cell(nPop*totalCases,1);
    activitiesPH2 = cell(nPop*totalCases,1);
    activitiesPH3 = cell(nPop*totalCases,1);
    errors = zeros(nPop*totalCases,1);
    
    runData.ssCurves{runGeneration} = cell(nPop,totalCases);
    runData.err{runGeneration} = zeros(nPop,totalCases);
    
    %% Run write Sx File and immediately a.out in parallel
    parfor i = 1:size(targetFolders,1)
        % Assign current folder to curFold to make reading later code
        % easier
        curFold = targetFolders{i};
        % Find the upper path of the folder for later use
        [uppath,~] = fileparts(curFold);
        
        % Find the difference in string size
        pathLenDiff = length(curFold) - length(uppath); % this is the length if /C1_1... /C1_100
        % Difference is used to properly index par, end-pathLenDiff+5
        % returns the numbers at the end of curFold: the 100 in /C1_100
        parInd = str2double(curFold(end-pathLenDiff+5:end));
        
        curInputs=inputs;
        curInputs(toRefine,1)= par(parInd,:)';
        
        % Move to optimization folder to run writeSxFile
%         cd([uppath,'/Optimization']); % Moved to main
        WriteSxFile([curFold,'/Titanium.sx'],curInputs,linCases{i});
        WriteLoadFile([curFold,'/load.pro'],systemParams,linLoading{i})
%         WriteLoadFile([curFold,'/Titanium.sx'],MaxStrain,)
        % Move into the EPSC folder to run a.out or the windows exe
        cd(targetFolders{i});
        if (isystem == 0)
            % run the compression files in batch for windows
            [~,~] = system('epsc4_6_4_vol_avg_spin.exe');
        elseif (isystem == 1)
            % run the compression files in batch for linux
            [~,~] = system([curFold,'/a.out']);
        end
        cd(rootPath);
        % Load in the freshly baked epsc3.out
        vpscData = importdata([curFold,'/STR_STR.OUT']);
        vpscData = vpscData.data;
        curSSCurves{i} = vpscData;
        
        % Calculate error ------------------------------------------------
        % Prepare the simulated model
        curDir = str2double(curFold(end-pathLenDiff+3));
        VPSCSTRSTROffset = 2; % accounting for the first two columns being mean stress/strain
        simModx = abs(vpscData(:,curDir+VPSCSTRSTROffset));
        simMody = abs(vpscData(:,curDir+6+VPSCSTRSTROffset));
        
        % Determine the case
        modelCountCumSum = cumsum(modelCaseDist);
        caseNum = ceil(i/nPop);
        if (caseNum <= modelCountCumSum(1))
            curCase = 1;
            curModel = caseNum;
        else
            for curCase = 2:length(modelCountCumSum)
                if (caseNum > modelCountCumSum(curCase-1) && caseNum <= modelCountCumSum(curCase))
                    curModel = caseNum - modelCountCumSum(curCase-1);
                    break;
                end
            end
        end
        curmodeledCurves = load([rootPath,'/',expModNames{curCase}]);
        curModeledCurves = curmodeledCurves.curModeledCurves;
        origModel = curModeledCurves{curModel};
        origModelx = origModel(:,1); % Select the experimental models to be fit
        origModely = origModel(:,2);
        
        % Fit simulated model data to a function
        [simModelx, simModely] = prepareCurveData(simModx,simMody);
        [fitresult_mod, ~] = fit(simModelx,simModely,'smoothingspline');
        
        % shift enabled implies it is not fitting a model, and
        % error comes from the shifted range
        if (ishift == 1)
            errors(i) = calcError(1,origModelx,origModely,fitresult_mod,fitRange,...
                10,-0.01,0.01);
        else % shift not enabled, could still be fitting a model or fitting part of the curve
            if (ifitpart == 0)
                % if fitting the whole curve
                errors(i) = calcError(2,origModelx,origModely,fitresult_mod,fitRange,simModelx);
            elseif (ifitpart == 1)
                % If only fitting a part of the curve
                errors(i) = calcError(3,origModelx,origModely,fitresult_mod,fitRange);
            end
        end
        
        vpscData = importdata([curFold,'/ACT_PH1.OUT']);
        vpscData = vpscData.data;
        activitiesPH1{i} = vpscData;
        vpscData = importdata([curFold,'/ACT_PH2.OUT']);
        vpscData = vpscData.data;
        activitiesPH2{i} = vpscData;
        vpscData = importdata([curFold,'/ACT_PH3.OUT']);
        vpscData = vpscData.data;
        activitiesPH3{i} = vpscData;
    end
    runData.ssCurves{runGeneration} = reshape(curSSCurves,[nPop,totalCases]);
    runData.activitiesPH1{runGeneration} = reshape(activitiesPH1,[nPop,totalCases]);
    runData.activitiesPH2{runGeneration} = reshape(activitiesPH2,[nPop,totalCases]);
    runData.activitiesPH3{runGeneration} = reshape(activitiesPH3,[nPop,totalCases]);
    runData.err{runGeneration} = reshape(errors,[nPop,totalCases]);
    errTotal = reshape(errors,[nPop,totalCases]);
    errTotal = sum(errTotal,2);
    % Delete all files and return back to the root directory where the
    % data is stored stored
    cd(rootPath);
end