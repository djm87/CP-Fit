function optimizationStart(cases,fitParam,info,GAoptions,runFoldName)
% Read in fitting parameters for DD law
Table = fitParam;
% Determine the fitting recipe
toRefine = find(Table.fitFlag == true);

systemParams = {cases;fitParam;info;GAoptions;runFoldName};

%% Global Variables
global runData;
runData.params = cell(GAoptions.MaxGenerations,1); % number of generations
runData.mainSimData = cell(GAoptions.MaxGenerations,1); % number of generations
% runData.activitiesPH1 = cell(GAoptions.MaxGenerations,1); % number of generations
% runData.activitiesPH2 = cell(GAoptions.MaxGenerations,1); % number of generations
% runData.activitiesPH3 = cell(GAoptions.MaxGenerations,1); % number of generations
runData.err = cell(GAoptions.MaxGenerations,1); % number of generations
global runGeneration;
runGeneration = 0;

%% Optimization Run
if (info.GAinp.iSingleRunGA == 0)
    fun = @modelWrapErrFunction;
    
    lb = Table{toRefine,{'lowerBound'}};
    ub = Table{toRefine,{'upperBound'}};     
    nVars = length(ub);
    A = info.GAinp.A;
    b = info.GAinp.b;
    Aeq = info.GAinp.Aeq;
    beq = info.GAinp.beq;
    nonlcon = info.GAinp.nonlcon;
    IntCon = info.GAinp.IntCon(1):info.GAinp.IntCon(2);
    
    if (info.GAinp.gaType == 1)
        [xOut,Fval,exitFlag,Output] = ga(@(x) fun(x,systemParams),...
            nVars,A,b,Aeq,beq,lb,ub,nonlcon,IntCon,GAoptions);
    else
        [xOut,Fval,exitFlag,Output] = gamultiobj(@(x) fun(x,systemParams),...
            nVars,A,b,Aeq,beq,lb,ub,nonlcon,GAoptions);
    end
    
%     FigList = findobj(allchild(0), 'flat', 'Type', 'figure');
%     savefig(FigList,'GAFigs.fig');
    
    if (info.fitStrat.isave == 1)
        %Addded bit so there is no overwriting
        fname1='AllVars.mat';
        fname2='runData.mat';
        fname3='Fitting_inputs.csv';
        cnt=0;
        
        while (exist(fname1,'file'))
            cnt=cnt+1;
            fname1=['AllVars' int2str(cnt) '.mat'];
            fname2=['runData' int2str(cnt) '.mat'];
            fname3=['Fitting_inputs' int2str(cnt) '.csv'];
        end
        save(fname1, '-regexp', '^(?!(runData)$).');
        save(fname2,'runData');
    end
else
    fname4 = 'oneRunAllVars.mat';
    fname5 = 'oneRunRunData.mat';
    errTotal = modelWrapErrFunction(Table{toRefine,{'Parameters'}}',systemParams);
    save(fname4,'-regexp', '^(?!(runData)$).');
    save(fname5,'runData');
end

%% Make plot of the model of the lowest error in the last generation
plotGen = GAoptions.MaxGenerations - 1;
[~,plotPop] = min(sum(runData.err{plotGen,1},2));

Table{toRefine,{'Parameters'}}=runData.params{plotGen,1}(plotPop,:)';
Table{toRefine,{'lowerBound'}}=lb;
Table{toRefine,{'upperBound'}}=ub;
Table{1:end,{'fitFlag'}}=0;
Table{toRefine,{'fitFlag'}}=1;
writetable(Table,fname3);

caseDataFiles = cases{1:end,'FilePath'};
if (info.fitStrat.iplot == 2)
    for i = 1:numel(caseDataFiles)
        exp = importdata(caseDataFiles{i});
        sim = runData.mainSimData{plotGen}{plotPop,i};
%         ActPH1=runData.activitiesPH1{plotGen}{plotPop,i};
%         ActPH2=runData.activitiesPH2{plotGen}{plotPop,i};
%         ActPH3=runData.activitiesPH3{plotGen}{plotPop,i};
        
        figure(cases.PlotFigure(i));
        subplot(1,2,1);
        hold on;
        plot(exp(:,1),exp(:,2),'--');
        plot(sim(:,1),exp(:,2));
        
%         strain=ActPH1(:,1);
%         phaseFrac=ActPH1(:,2);
%         subplot(1,2,2);
%         hold on;
%         plot(strain,smooth(ActPH1(:,4),5),'-','LineWidth',3, 'Color','b')
%         hold on
%         plot(strain,smooth(ActPH1(:,6),5),'-','LineWidth',3, 'Color','g')
%         plot(strain,smooth(ActPH1(:,5),5),'-','LineWidth',3, 'Color','r')
%         plot(strain,phaseFrac,'--','lineWidth',3, 'Color','k')
%         plot(strain,smooth(ActPH1(:,9),5),'--','lineWidth',3, 'Color','y')
%         plot(strain,smooth(ActPH1(:,10),5),'--','lineWidth',3, 'Color','c')
    end
end
end