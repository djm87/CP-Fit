function optimizationStart(cases,fitParam,info,runSource,GAoptions)

% Read in fitting parameters for DD law
Table = fitParam;
% Determine the fitting recipe
toRefine = find(Table.fitFlag == true);

%% Global Variables
global runData;
runData.params = cell(GAoptions.MaxGenerations,1); % number of generations
runData.mainSimData = cell(GAoptions.MaxGenerations,1); % number of generations
runData.activitiesPH1 = cell(GAoptions.MaxGenerations,1); % number of generations
runData.activitiesPH2 = cell(GAoptions.MaxGenerations,1); % number of generations
runData.activitiesPH3 = cell(GAoptions.MaxGenerations,1); % number of generations
runData.err = cell(GAoptions.MaxGenerations,1); % number of generations
global runGeneration;
runGeneration = 0;

%% Optimization Run
if info.GAinp.iSingleRunGA==true
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
        [xOut,Fval,exitFlag,Output] = ga(@(x) fun(x,cases,fitParam,info),...
            nVars,A,b,Aeq,beq,lb,ub,nonlcon,IntCon,GAoptions);
    else
        [xOut,Fval,exitFlag,Output] = gamultiobj(@(x) fun(x,cases,fitParam,info),...
            nVars,A,b,Aeq,beq,lb,ub,nonlcon,GAoptions);
    end
    
    if (info.fitStrat.isave == 1)
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
    end
else 
    errTotal = VPSC_WrapperFunction_MultObj(Table{toRefine,{'Parameters'}}',systemParams);
end

%% Make plot of the model 
% if runFit==true
%     gen=8;
%     [~,idx]=min(sum(runData.err{gen,1},2))
% %     [~,idx]=min(runData.err{gen,1}(:,7))
%     Table{toRefine,{'Parameters'}}=runData.params{gen,1}(idx,:)';
%     Table{toRefine,{'lowerBound'}}=lb;
%     Table{toRefine,{'upperBound'}}=ub;
% 
% 
%     Table{1:end,{'fitFlag'}}=0;
%     Table{toRefine,{'fitFlag'}}=1;
% 
%     writetable(Table,fname3);
% else
%    gen=1;
%     idx=1; 
% end
%    
% 
% shift_SS=[-0.0,-0.,-0.0,0.0,0.0,0.0,0.0,0.0,0.0]
% for i = [gen]%1:runGeneration %Generation loops
%     for j = [idx]%1:size(runData.activitiesPH1{1,1},1) %population loop
%         for k= 1:length(dataFileNames)
%             
%         SSm=runData.ssCurves{i,1}{j,k};
%         SSx=load([rootPath,'/YieldSearch/',dataFileNames{k}]);
%         ActPH1=runData.activitiesPH1{i,1}{j,k};
%         ActPH2=runData.activitiesPH2{i,1}{j,k};
%         ActPH3=runData.activitiesPH3{i,1}{j,k};
%         
%         figure('Name',['case: ' int2str(j) ', data name: ' dataFileNames{k}],...
%             'Position',[188 541 1000 420])
%         subplot(1,2,1)
%         plot(SSx(:,1)+shift_SS(k),SSx(:,2),'LineStyle','-','LineWidth',2, 'Color','b')
%         hold on
%         plot(SSm(:,1),SSm(:,2),'LineStyle','--','LineWidth',2, 'Color','r')
% 
% 
%         xlabel('True Strain ','FontSize',14)
%         ylabel('True Stress [MPa] ','FontSize',14)
%         % title(['$$ \rm Compression $$'] ,'FontSize',20,'interpreter','latex')
%         % set(gca,'XTick',[0 0.025 0.05 0.075 0.1 0.125 0.15 0.175 0.2 0.225 0.25 0.275 0.3 0.325 0.35 0.375 0.4],'FontSize',16)
%         % set(gca,'XTick',[0 0.05 0.1 0.15 0.2 0.25],'FontSize',18)
%         % set(gca,'YTick',[100 200 300 400 500],'FontSize',18)
%         
% %          leg1=legend('Exp.','Model');%)%,'TT1')%,'TT1 VF')%'Pyramidal-1st order',
% %         set(leg1,'Orientation','vertical','Location','best',...
% %             'FontSize',12);
%         axis tight
%         box on
%         grid off
%         xlim([0 fitRange(2)]) 
%         ylim([0,900])
%         
%         subplot(1,2,2)
%         strain=ActPH1(:,1);
%         phaseFrac=ActPH1(:,2);
%         plot(strain,smooth(ActPH1(:,4),5),'-','LineWidth',3, 'Color','b')
%         hold on
%         plot(strain,smooth(ActPH1(:,6),5),'-','LineWidth',3, 'Color','g')
%         plot(strain,smooth(ActPH1(:,5),5),'-','LineWidth',3, 'Color','r')
%         plot(strain,phaseFrac,'--','lineWidth',3, 'Color','k')
%         plot(strain,smooth(ActPH1(:,9),5),'--','lineWidth',3, 'Color','y')
%         plot(strain,smooth(ActPH1(:,10),5),'--','lineWidth',3, 'Color','c')
% 
%         xlim([0 fitRange(2)]) 
%         ylim([0 1])
%         yticks([0,0.5,1])
%         yticklabels({'0','0.5','1'})
% %         xticks([0,0.3,0.6,0.9])
% %         xticklabels({'0','0.3','0.6','0.9'})
%         title('Parent Activities','FontSize',14) 
%         xlabel('True Strain','FontSize',14) 
% % %                 leg1=legend('<a> prism','<a> basal','<c+a> pyr I','Phase Frac. Act','T VF','C VF');%)%,'TT1')%,'TT1 VF')%'Pyramidal-1st order',
% %         set(leg1,'Orientation','vertical','Location','best',...
% %             'FontSize',12);
%         ylabel('Activity','FontSize',14) 
%         % set(gca,'XTick',[0 0.05 0.1 0.15 0.2 0.25],'FontSize',18)
%         % set(gca,'YTick',[0 0.25 0.5 0.75 1],'FontSize',12)
%         % axis tight
%         % box off
%         % grid off
% 
% %         leg1=legend('<a> prism','<a> basal','<c+a> pyr II','Par. VF');%)%,'TT1')%,'TT1 VF')%'Pyramidal-1st order',
% %         set(leg1,'Orientation','horizontal','Location','southoutside',...
% %             'FontSize',12);
% %         legend('boxoff')
% %         
% %        subplot(1,3,3)
% %         strain=ActPH2(:,1);
% %         phaseFrac=ActPH2(:,2);
% %         plot(strain,smooth(ActPH2(:,4),5),'-','LineWidth',3, 'Color','b')
% %         hold on
% %         plot(strain,smooth(ActPH2(:,6),5),'-','LineWidth',3, 'Color','g')
% %         plot(strain,smooth(ActPH2(:,5),5),'-','LineWidth',3, 'Color','r')
% %         plot(strain,phaseFrac,'--','lineWidth',3, 'Color','k')
% %         plot(strain,smooth(ActPH2(:,5),9),'--','lineWidth',3, 'Color','y')
% %         plot(strain,smooth(ActPH2(:,5),10),'--','lineWidth',3, 'Color','c')
% % 
% %         xlim([0 0.4]) 
% %         ylim([0 1])
% %         yticks([0,0.5,1])
% %         yticklabels({'0','0.5','1'})
% % %         xticks([0,0.3,0.6,0.9])
% % %         xticklabels({'0','0.3','0.6','0.9'})
% %         title('Twin Activities') 
% %         xlabel('True Strain')
% %         ylabel('Activity')
% %         % set(gca,'XTick',[0 0.05 0.1 0.15 0.2 0.25],'FontSize',18)
% %         % set(gca,'YTick',[0 0.25 0.5 0.75 1],'FontSize',12)
% %         % axis tight
% %         % box off
% %         % grid off
% % 
% %         leg1=legend('<a> prism','<a> basal','<c+a> pyr II','Phase F Act','T VF','C VF');%)%,'TT1')%,'TT1 VF')%'Pyramidal-1st order',
% %         set(leg1,'Orientation','vertical','Location','best',...
% %             'FontSize',12);
% %         legend('boxoff') 
%         print( '-opengl', '-dtiff', '-r300',['Result_' int2str(k)]); %saveas(h,'odf_ND_ini','pdf');
% 
%         end
%     end
% end

end