% global runData;
% global runGeneration;

%% Plot the ratios of the three CRSS
% data = xOut(:,1:3);
% data = data./(data(:,1));
% xL = 1:size(data,1);
% figure;
% plot(xL,data);
% legend('prismatic','basal','pyramidal');

%% Prepare new folders with results


%% Plot the simulated results
% You can run this section directly after the main process is completed. Do
% not clear variables or figures.
load('AllVars.mat','xOut');
load('AllVars.mat','cases');
load('AllVars.mat','dataFileNames');
load('AllVars.mat','rootPath');
xOpt = xOut;
% Write new titanium files using xOut values
rtPath = rootPath;
test12Path = 'E:/Ti/Results/Result10/12um';
test20Path = 'E:/Ti/Results/Result10/20um';
test30Path = 'E:/Ti/Results/Result10/30um';

allPaths = {[test12Path,'/',cases{1}(1,:),'_1'];...
    [test12Path,'/',cases{1}(2,:),'_1'];...
    [test12Path,'/',cases{1}(3,:),'_1'];...
            [test20Path,'/',cases{2}(1,:),'_1'];...
    [test20Path,'/',cases{2}(2,:),'_1'];
    [test20Path,'/',cases{2}(3,:),'_1'];
    [test20Path,'/',cases{2}(4,:),'_1'];
    [test20Path,'/',cases{2}(5,:),'_1'];
    [test20Path,'/',cases{2}(6,:),'_1'];
            [test30Path,'/',cases{3}(1,:),'_1'];
    [test30Path,'/',cases{3}(2,:),'_1']};

err = zeros(size(allPaths,1),1);
ft = fittype( 'smoothingspline' );

for i = 1:size(xOpt,1)
    curX = xOpt(i,:);
    parfor j = 1:size(allPaths,1)
        cd(allPaths{j});
        WriteSxFile('Titanium.sx',curX);
        [~,~] = system('epsc4_6_4_vol_avg_spin.exe');
    end
    cd('E:/Ti/Results/Result10');
    for j = 1:size(allPaths,1)
        mod = importdata([allPaths{j},'/epsc3.out']);
        curMod = mod.data;
        curDir = str2double(allPaths{j}(end-2));
        
        modx = abs(curMod(:,curDir));
        mody = abs(curMod(:,curDir+6));
        
        orgMod = importdata(dataFileNames{j});
        orgx = abs(orgMod(:,1));
        orgy = abs(orgMod(:,2));
        
        figure(j);
        hold on;
        plot(modx,mody,'b--');
        plot(orgx,orgy,'k');
        title(['Optimized Model and Original Model ',dataFileNames{j}(4:5),' ',dataFileNames{j}(7:8)]);
        xlabel('Strain');
        ylabel('Stress (MPa)');
        legend('Location','best','optimized','original');
        xlim([0 0.05]);
        ylim([0 400]);
        saveas(gcf,['Result10Zoomed_',int2str(j),'.png']);
        
        % error analysis
        [xData, yData] = prepareCurveData(modx,mody);
        [fitresult_mod, ~] = fit(xData,yData,ft);
        
        [xData_o, yData_o] = prepareCurveData(orgx,orgy);
        [fit_o, ~] = fit(xData_o,yData_o,ft);
        
        xIn = linspace(0,0.04,100);
        
        funDiff = fitresult_mod(xIn) - fit_o(xIn);
        diffSq = funDiff.^2;
        diffSqSum = sum(diffSq);
        err(j) = sqrt(diffSqSum);
    end
end

% cd('/media/data/Jesse/Ti/Model/InitialYield_New/TestResult3');
% count = 1;
% cCount = 1;
% for i = 1:11
%     figure(i);
%     if (count == 4 && cCount == 1)
%         cCount = cCount + 1;
%         count = 1;
%     elseif (count == 7 && cCount == 2)
%         cCount = cCount + 1;
%         count = 1;
%     end
%     saveas(gcf,['Result_',models{cCount},'_',cases{cCount}(count,:),'.png']);
%     count = count + 1;
% end