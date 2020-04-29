function optimizationStart(cases,fitParam,fitRecipe,info,GAoptions,runFoldName,PFFits,expData)
systemParams = {cases;fitParam;fitRecipe;info;GAoptions;runFoldName;PFFits;expData};

%% Global Variables
global runData;
runData.params = cell(GAoptions.MaxGenerations,1); % number of generations
runData.lowestErrSimData = cell(GAoptions.MaxGenerations,1);
runData.shiftinds = cell(GAoptions.MaxGenerations,1);
runData.err = cell(GAoptions.MaxGenerations,1);
global runGeneration;
runGeneration = 0;

%% Optimization Run
if (info.GAinp.iSingleRunGA == 0)
    fun = @modelWrapErrFunction;
    
    lb = fitParam{fitRecipe,{'lowerBound'}};
    ub = fitParam{fitRecipe,{'upperBound'}};     
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
    errTotal = modelWrapErrFunction(fitParam{fitRecipe,{'Parameters'}}',systemParams);
    save(fname4,'-regexp', '^(?!(runData)$).');
    save(fname5,'runData');
end

%% Make plot of the model of the lowest error in the last generation
[~,plotPop] = min(sum(runData.err{runGeneration,1},2));

fitParam{fitRecipe,{'Parameters'}}=runData.params{runGeneration,1}(plotPop,:)';
writetable(fitParam,fname3);

caseDataFiles = cases{1:end,'FilePath'};
if (info.fitStrat.iplot == 1)
    figure;
    shift = runData.shiftinds{runGeneration};
    shiftAmt = linspace(info.fitStrat.ishift(2),info.fitStrat.ishift(3),info.fitStrat.ishift(4));
    for i = 1:numel(caseDataFiles)
        exp = importdata(caseDataFiles{i});
        sim = runData.lowestErrSimData{runGeneration}{1,i};
        
        subplot(ceil(length(caseDataFiles)/3),3,cases.PlotFigure(i));
        
        hold on;
        plot(exp(:,1),exp(:,2),'r','LineWidth',2);
        plot(sim(:,1)+shiftAmt(shift(i)),sim(:,2),'b--','LineWidth',2);
        xlim([-0.2 0]);
        xlabel('True strain');
        ylabel('True stress (MPa)');
    end
elseif (info.fitStrat.iplot == 2)
    for i = 1:numel(caseDataFiles)
        exp = importdata(caseDataFiles{i});
        sim = runData.lowestErrSimData{runGeneration}{1,i};
%         ActPH1=runData.activitiesPH1{plotGen}{plotPop,i};
%         ActPH2=runData.activitiesPH2{plotGen}{plotPop,i};
%         ActPH3=runData.activitiesPH3{plotGen}{plotPop,i};
        figure(cases.PlotFigure(i));
        hold on;
        plot(exp(:,1),exp(:,2),'r--','LineWidth',2);
        plot(sim(:,1),sim(:,2),'b-','LineWidth',2);
        xlabel('True strain');
        ylabel('True stress (MPa)');
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