cd 'E:\3 Ti\Outputs\1 Ti Outputs\2 Set2\3 Neff Fitting 2 points shift reduced ranges'
%%
cd 'E:\3 Ti\CP-Fit\GAFit'
caseIDs = cases{1:end,'CaseIdentifier'};
lowestErr = 10000;
lowestErrPop = 0;
lowestErrGen = 0;
for i = 2:runGeneration
    [curMinErr,plotPop] = min(runData.err{i,1}(:,5));
    if (curMinErr < lowestErr)
        lowestErr = curMinErr;
        lowestErrPop = plotPop;
        lowestErrGen = i;
    end
end
caseDataFiles = cases{1:end,'FilePath'};
[~,plotPop1] = min(sum(runData.err{runGeneration,1}(:,11:end),2))
[~,plotPop2] = min(sum(runData.err{runGeneration,1},2))
figure;
shift = runData.shiftinds{runGeneration};
shiftAmt = linspace(info.fitStrat.ishift(2),info.fitStrat.ishift(3),info.fitStrat.ishift(4));
for i = 1:numel(caseDataFiles)
    expData2 = importVPSCout(caseDataFiles{i},0);
    sim = runData.lowestErrSimData{runGeneration}{1,i};%importVPSCout(['E:\3 Ti\Outputs\1 Ti Outputs\2 Set2\brad_625_800 10 cases recipe 2 from try 2\265\',num2str(i),'\STR_STR.OUT'],1);
    
    subplot(ceil(length(caseDataFiles)/3),3,cases.PlotFigure(i));
    
    hold on;
    plot(expData2(:,1),expData2(:,2),'r','LineWidth',2);
    plot(sim(:,1),sim(:,2),'b--','LineWidth',2);%+shiftAmt(shift(i))
    xlim([-0.4 0]);
    xlabel('True strain');
    ylabel('True stress (MPa)');
end

figure;
for i = 1:numel(caseDataFiles)
    sim = runData.lowestErrVFData{runGeneration}{1,i};%importVPSCout(['E:\3 Ti\Outputs\1 Ti Outputs\2 Set2\3 Neff Fitting 2 points shift reduced ranges\11 results like 8 but with different ranges and fitting texture\193\',num2str(i),'\ACT_PH1.OUT'],1);
    
    subplot(ceil(length(caseDataFiles)/3),3,cases.PlotFigure(i));
    
    hold on;
    plot(sim(:,1),sim(:,2),'k','LineWidth',2);
    plot(sim(:,1),sim(:,3),'k-.','LineWidth',2);
    plot(sim(:,1),sim(:,4),'k--','LineWidth',2);
    plot(sim(:,1),sim(:,5),'b--','LineWidth',2);
    plot(sim(:,1),sim(:,6),'r','LineWidth',2);
%     xlim([-0.2 0]);
    xlabel('True strain');
    ylabel('Activities and VF');
    box off;
end
legend('Pri','Pyr','Bas','TT','CT','Location','EastOutside');
legend('box','off');

VFExp = {[3.08,54.96,1.2;0.07,1.03,5.14],...
    [0.27,6.76,0.72;2.53,7.53,0.35],...
    [0.15,17.31,0.82;0.36,7.42,0.18],...
    [6.59,46.34,3.18;0.29,12.81,13.43],...
    [0.46;2.38],...
    [1.07,10.71,0.05,9.93;9.98,33.46,0.02,3.74]};
figure;
for i = 1:numel(caseDataFiles)
    sim = runData.lowestErrVFData{runGeneration}{1,i};%importVPSCout(['E:\3 Ti\Outputs\1 Ti Outputs\2 Set2\3 Neff Fitting 2 points shift reduced ranges\11 results like 8 but with different ranges and fitting texture\193\',num2str(i),'\ACT_PH1.OUT'],1);
    sim2 = runData.lowestErrVFData2{runGeneration}{1,i};%importVPSCout(['E:\3 Ti\Outputs\1 Ti Outputs\2 Set2\3 Neff Fitting 2 points shift reduced ranges\11 results like 8 but with different ranges and fitting texture\193\',num2str(i),'\ACT_PH1.OUT'],1);
    subplot(ceil(length(caseDataFiles)/3),3,cases.PlotFigure(i));
    
    hold on;
    plot(sim(:,1),sim(:,5),'k','LineWidth',2);
    plot(sim(:,1),sim(:,6),'b','LineWidth',2);
    plot(sim2(:,1),sim2(:,2),'r','LineWidth',2);
    plot(sim2(:,1),sim2(:,3),'m','LineWidth',2);
    if (i >= 5 && i < 9)
        plot(0.05,VFExp{i-4}(1,1)/100,'kO');
        plot(0.2,VFExp{i-4}(1,2)/100,'kO');
        plot(0.2,VFExp{i-4}(1,3)/100,'rO');
        plot(0.05,VFExp{i-4}(2,1)/100,'bO');
        plot(0.2,VFExp{i-4}(2,2)/100,'bO');
        plot(0.2,VFExp{i-4}(2,3)/100,'mO');
    elseif (i == 9)
        plot(0.05,VFExp{i-4}(1,1)/100,'kO');
        plot(0.05,VFExp{i-4}(2,1)/100,'bO');
    elseif (i == 10)
        plot(0.05,VFExp{i-4}(1,1)/100,'kO');
        plot(0.2,VFExp{i-4}(1,2)/100,'kO');
        plot(0.05,VFExp{i-4}(1,3)/100,'rO');
        plot(0.2,VFExp{i-4}(1,4)/100,'rO');
        plot(0.05,VFExp{i-4}(2,1)/100,'bO');
        plot(0.2,VFExp{i-4}(2,2)/100,'bO');
        plot(0.05,VFExp{i-4}(2,3)/100,'mO');
        plot(0.2,VFExp{i-4}(2,4)/100,'mO');
    end
%     xlim([-0.2 0]);
    xlabel('True strain');
    ylabel('Activities and VF');
    box off;
end
legend('TT','CT','TT2','CT2','Location','EastOutside');
legend('box','off');
%%
% scaling = fitParam.scaling(fitRecipe);
% sumErr = sum(runData.err{13,1},2);
% [a,ind] = sort(sumErr);
% t = zeros(length(b),5);
% b = [2.95,2.95,5.54,0.30174,0.2725]*1e-10;
% mu = 44e9;
% d = 25.5e-6;
% 
% for i = 1:length(ind)
%     curParams = runData.params{13,1}(ind(i),:);
%     
%     t0 = curParams([1 3 5 7 9]).*scaling([1 3 5 7 9])'; %MPa
%     HP =curParams([2 4 6 8 10]).*scaling([2 4 6 8 10])';
%     t(i,:) = t0 + mu.*HP.*sqrt(b./d)./1e6;
% end