%% Initilialize Matlab Environement 
% clear all;
close all;
clc;
format shortG
tic
%% Read in fitting parameters for DD law
Table=readtable('Fitting_inputs10_custom.csv');
inputs = Table{1:end,{'Parameters','lowerBound','upperBound','fitFlag','scaling'}};
find(Table{1:end,{'fitFlag'}}==true)
%% Select a Refinement Recipe 

recipe = 2;

% 1 Refine Tau0 and HP for slip 
if recipe==1; toRefine = [6,7,15,16,24,25];end

% 2 Refine Tau0 and HP for slip and twin 
if recipe==2;toRefine = [6,7,15,16,24,25,33,34,38,39];end

% 3 Refine Tau0 and HP for slip with k1
if recipe==3;toRefine = [6,7,8,15,16,17,24,25,26,33,34,38,39];end

% 4 Refine Tau0 and HP for slip and twin with k1 and g 
if recipe==4;toRefine = [6,7,8,9,15,16,17,18,24,25,26,27,33,34,38,39];end

% 4 Refine Tau0 and HP for slip with k1 and g 
if recipe==5;toRefine = [6,7,8,9,15,16,17,18,24,25,26,27];end

% 6 Refine DD law + twin HP/CRSS
if recipe==6;toRefine = [8,9,10,11,12,17,18,19,20,21,26,27,28,29,30];end

% 7 Refine select DD law + twin HP/CRSS
if recipe==7;toRefine = [8,9,11,17,18,20,26,27,29,33,34,38,39];end

% 8 Refine HP for twin with k1 and g
if recipe==8;toRefine = [8,9,17,18,26,27,33,34,38,39];end

%% Set environment and runtime options 
%GA parametrs - should be chosen carefully...
maxGens=9;
maxPop=10*length(toRefine);

runFit=true; %when false it runs a single case

boundsType=2; %1 for blanket percent difference, 2 for custom read in from csv 
boundDiff=0.6; %used if boundsType=1 

% Additional Options
ishift = 0; % 1 - enable shifting - if fitting a model this option is ignored
ifitpart = 1; % DS: no longer relevant 1 - calculates error on a small portion of the entire curve, 0 - calculates error on the whole curve (old
isystem = 1; % 0 - windows system running .exe, 1 - linux system running ./a.out
fitRange = [0.03,0.3]; %determines the level of strain in the model


%% Global Variables
global runData;
runData.params = cell(maxGens,1); % number of generations
runData.ssCurves = cell(maxGens,1); % number of generations
runData.activitiesPH1 = cell(maxGens,1); % number of generations
runData.activitiesPH2 = cell(maxGens,1); % number of generations
runData.activitiesPH3 = cell(maxGens,1); % number of generations
runData.err = cell(maxGens,1); % number of generations
global runGeneration; 
runGeneration = 0; 
%% Initialize the fitting cases

%DS:This section is very redundant..

% stress strain curve names
dataFileNames = {'exp12_C1.SS';'exp12_C2.SS';...%;'exp12_C3.SS';...
               'exp400_C1.SS';'exp400_C2.SS';'exp400_C3.SS';...
               'exp625_C1.SS';'exp625_C2.SS';'exp625_C3.SS'}; % Experimental SS curves

% The names of the .mat files that contains the experimental or modeled
% stress strain curves to be calculated for error
expModNames = {'exp_12um.mat';'exp_400C.mat';'exp_625C.mat'}; % matches the .mat files created by yieldSearch           
           
% Grain size cases or folders that are kept seperate. These Are folders in
% the main directory
models = {'12um';'400C';'625C'}; % Model folder names

% Cases that occur in each model 
cases = {['C1';'C2'];...
         ['C1';'C2';'C3'];...
         ['C1';'C2';'C3']}; % must match size of ss curves

% Creates the absolute paths for the models
rootPath = pwd;
modelPaths = cell(size(models));
for i = 1:size(models,1)
    modelPaths{i} = [rootPath,'/',models{i}];
end
     
modelCaseDist = zeros(size(cases));
for i = 1:size(cases,1)
    modelCaseDist(i) = size(cases{i},1); % makes parallel work
end

% Define the list of parameters to be passed into each function that needs
% the information.
systemParams = {models;...                  % systemParams{1}, list of models
                rootPath;...                % systemParams{2}, the root directory
                modelPaths;...              % systemParams{3}, name of the models
                cases;...                   % systemParams{4}, name of the cases
                0;...                       % systemParams{5}, used by setup folders to know how many to setup (dummy input)
                modelCaseDist;...           % systemParams{6}, how many cases in each model
                expModNames;...             % systemParams{7}, the names of the files containing the experimental or modeled data
                fitRange;...                % systemParams{8}, the shifting plan for the experimental data
                [ishift,ifitpart,isystem];... % systemParams{9}, the additional options
                toRefine;...                % systemParams{10}, which parameters are being changed
                inputs};                     % systemParams{11}, All the Input information
%% Initialize stress strain data

% DS:Does this do anything now?

% Call the yieldSearch function to initialize the stress strain data
cd([rootPath,'/YieldSearch']);
yieldSearch(dataFileNames,systemParams);
cd(rootPath);
%% Optimization Run
% fun = @EPSC_WrapperFunction_MultObj;
if runFit==true
    fun = @VPSC_WrapperFunction_MultObj;
    % 
    if boundsType==1
        lb = (1-boundDiff).*Table{toRefine,{'Parameters'}};
        ub = (1+boundDiff).*Table{toRefine,{'Parameters'}};  
    else 
        lb = Table{toRefine,{'lowerBound'}};
        ub = Table{toRefine,{'upperBound'}};     
    end
    nVars = length(ub);
    % ub = [100,28r0,280,2.5,2.5,2.5];%,30]; % grain sizes commented out
    A = [];
    b = [];
    Aeq = [];
    beq = [];

    % For integer parameters replace the following until the next ------------
    options = optimoptions('gamultiobj',...
        'PlotFcn',{@gaplotrange,@gaplotselection,@gaplotbestindiv...
        @gaplotscorediversity,@gaplotstopping,@gaplotgenealogy},...
        'FunctionTolerance',0.1,...
        'MaxGenerations',maxGens,...
        'MaxStallGenerations',20,...
        'PopulationSize',maxPop,...
        'ParetoFraction',0.35,...
        'Display','iter',...
        'UseVectorized',true);
    options.DistanceMeasureFcn = {@distancecrowding,'genotype'};

    [xOut,Fval,exitFlag,Output] = gamultiobj(@(x) fun(x,systemParams),nVars,A, ...
        b,Aeq,beq,lb,ub,options);
    % ------------------------------------------------------------------------
    % [x,fval,exitflag] = ga(@(x) fun(x,systemParams),nvars,A,b,[],[],...
    %     lb,ub,[],[1:6],options);
    % ------------------------------------------------------------------------

    %Addded bit so there is no overwriting
    cd(rootPath);
    fname1='AllVars.mat';
    fname2='runData.mat';
    fname3='Fitting_inputs.csv';
    cnt=0;

    while exist(fname1)~=false
        cnt=cnt+1;
        fname1=['AllVars' int2str(cnt) '.mat'];
        fname2=['runData' int2str(cnt) '.mat'];
        fname3=['Fitting_inputs' int2str(cnt) '.csv'];

    end
    save(fname1, '-regexp', '^(?!(runData)$).');
    save([rootPath,'/',fname2],'runData');
else 
    errTotal = VPSC_WrapperFunction_MultObj(Table{toRefine,{'Parameters'}}',systemParams);
end
%% Delete all the folders to save up space
% Read in the .mat file that contains the array of folders
% This assumes there are no subfolders since MATLAB cannot delete a folder
% with contents inside and there should be no folders in a normal operation
targetFolders = load('del_fold.mat');
targetFolders = targetFolders.targetFolders;
parfor i = 1:size(targetFolders,1)
    cd(targetFolders{i});
    
    dinfo = dir(pwd);
    for j = 3:size(dinfo,1)
        if (dinfo(j).isdir ~= 1)
            filename = fullfile(pwd,dinfo(j).name);
            delete(filename);
        end
    end
end
cd(rootPath);
for i=1:size(targetFolders,1)
    rmdir(targetFolders{i}); % remove all the now-empty folders
end

%% Make plot of the model 
if runFit==true
    gen=8;
    [~,idx]=min(sum(runData.err{gen,1},2))
%     [~,idx]=min(runData.err{gen,1}(:,7))
    Table{toRefine,{'Parameters'}}=runData.params{gen,1}(idx,:)';
    Table{toRefine,{'lowerBound'}}=lb;
    Table{toRefine,{'upperBound'}}=ub;


    Table{1:end,{'fitFlag'}}=0;
    Table{toRefine,{'fitFlag'}}=1;

    writetable(Table,fname3);
else
   gen=1;
    idx=1; 
end
   

shift_SS=[-0.0,-0.,-0.0,0.0,0.0,0.0,0.0,0.0,0.0]
for i = [gen]%1:runGeneration %Generation loops
    for j = [idx]%1:size(runData.activitiesPH1{1,1},1) %population loop
        for k= 1:length(dataFileNames)
            
        SSm=runData.ssCurves{i,1}{j,k};
        SSx=load([rootPath,'/YieldSearch/',dataFileNames{k}]);
        ActPH1=runData.activitiesPH1{i,1}{j,k};
        ActPH2=runData.activitiesPH2{i,1}{j,k};
        ActPH3=runData.activitiesPH3{i,1}{j,k};
        
        figure('Name',['case: ' int2str(j) ', data name: ' dataFileNames{k}],...
            'Position',[188 541 1000 420])
        subplot(1,2,1)
        plot(SSx(:,1)+shift_SS(k),SSx(:,2),'LineStyle','-','LineWidth',2, 'Color','b')
        hold on
        plot(SSm(:,1),SSm(:,2),'LineStyle','--','LineWidth',2, 'Color','r')


        xlabel('True Strain ','FontSize',14)
        ylabel('True Stress [MPa] ','FontSize',14)
        % title(['$$ \rm Compression $$'] ,'FontSize',20,'interpreter','latex')
        % set(gca,'XTick',[0 0.025 0.05 0.075 0.1 0.125 0.15 0.175 0.2 0.225 0.25 0.275 0.3 0.325 0.35 0.375 0.4],'FontSize',16)
        % set(gca,'XTick',[0 0.05 0.1 0.15 0.2 0.25],'FontSize',18)
        % set(gca,'YTick',[100 200 300 400 500],'FontSize',18)
        
%          leg1=legend('Exp.','Model');%)%,'TT1')%,'TT1 VF')%'Pyramidal-1st order',
%         set(leg1,'Orientation','vertical','Location','best',...
%             'FontSize',12);
        axis tight
        box on
        grid off
        xlim([0 fitRange(2)]) 
        ylim([0,900])
        
        subplot(1,2,2)
        strain=ActPH1(:,1);
        phaseFrac=ActPH1(:,2);
        plot(strain,smooth(ActPH1(:,4),5),'-','LineWidth',3, 'Color','b')
        hold on
        plot(strain,smooth(ActPH1(:,6),5),'-','LineWidth',3, 'Color','g')
        plot(strain,smooth(ActPH1(:,5),5),'-','LineWidth',3, 'Color','r')
        plot(strain,phaseFrac,'--','lineWidth',3, 'Color','k')
        plot(strain,smooth(ActPH1(:,9),5),'--','lineWidth',3, 'Color','y')
        plot(strain,smooth(ActPH1(:,10),5),'--','lineWidth',3, 'Color','c')

        xlim([0 fitRange(2)]) 
        ylim([0 1])
        yticks([0,0.5,1])
        yticklabels({'0','0.5','1'})
%         xticks([0,0.3,0.6,0.9])
%         xticklabels({'0','0.3','0.6','0.9'})
        title('Parent Activities','FontSize',14) 
        xlabel('True Strain','FontSize',14) 
% %                 leg1=legend('<a> prism','<a> basal','<c+a> pyr I','Phase Frac. Act','T VF','C VF');%)%,'TT1')%,'TT1 VF')%'Pyramidal-1st order',
%         set(leg1,'Orientation','vertical','Location','best',...
%             'FontSize',12);
        ylabel('Activity','FontSize',14) 
        % set(gca,'XTick',[0 0.05 0.1 0.15 0.2 0.25],'FontSize',18)
        % set(gca,'YTick',[0 0.25 0.5 0.75 1],'FontSize',12)
        % axis tight
        % box off
        % grid off

%         leg1=legend('<a> prism','<a> basal','<c+a> pyr II','Par. VF');%)%,'TT1')%,'TT1 VF')%'Pyramidal-1st order',
%         set(leg1,'Orientation','horizontal','Location','southoutside',...
%             'FontSize',12);
%         legend('boxoff')
%         
%        subplot(1,3,3)
%         strain=ActPH2(:,1);
%         phaseFrac=ActPH2(:,2);
%         plot(strain,smooth(ActPH2(:,4),5),'-','LineWidth',3, 'Color','b')
%         hold on
%         plot(strain,smooth(ActPH2(:,6),5),'-','LineWidth',3, 'Color','g')
%         plot(strain,smooth(ActPH2(:,5),5),'-','LineWidth',3, 'Color','r')
%         plot(strain,phaseFrac,'--','lineWidth',3, 'Color','k')
%         plot(strain,smooth(ActPH2(:,5),9),'--','lineWidth',3, 'Color','y')
%         plot(strain,smooth(ActPH2(:,5),10),'--','lineWidth',3, 'Color','c')
% 
%         xlim([0 0.4]) 
%         ylim([0 1])
%         yticks([0,0.5,1])
%         yticklabels({'0','0.5','1'})
% %         xticks([0,0.3,0.6,0.9])
% %         xticklabels({'0','0.3','0.6','0.9'})
%         title('Twin Activities') 
%         xlabel('True Strain')
%         ylabel('Activity')
%         % set(gca,'XTick',[0 0.05 0.1 0.15 0.2 0.25],'FontSize',18)
%         % set(gca,'YTick',[0 0.25 0.5 0.75 1],'FontSize',12)
%         % axis tight
%         % box off
%         % grid off
% 
%         leg1=legend('<a> prism','<a> basal','<c+a> pyr II','Phase F Act','T VF','C VF');%)%,'TT1')%,'TT1 VF')%'Pyramidal-1st order',
%         set(leg1,'Orientation','vertical','Location','best',...
%             'FontSize',12);
%         legend('boxoff') 
        print( '-opengl', '-dtiff', '-r300',['Result_' int2str(k)]); %saveas(h,'odf_ND_ini','pdf');

        end
    end
end
toc