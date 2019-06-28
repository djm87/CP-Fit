% exp yield at about 242 at 0.05362 and 84 index
% mod reach ~242 at 0.0067 and 68 index
exp12_C1 = importdata('exp12_C1.SS');
exp12_C1x = exp12_C1(:,1);
exp12_C1y = exp12_C1(:,2);
mod12_C1 = importdata('mod12_C1.SS');
mod12_C1x = mod12_C1(:,1);
mod12_C1y = mod12_C1(:,2);

figure(1);
hold on;
plot(exp12_C1x,exp12_C1y,'--');
plot(mod12_C1x,mod12_C1y);
xlim([0,0.02])
%
% C1 12um exp data edit
% slope = mod12_C1(10,2)/mod12_C1(10,1);
% yieldStrain = exp12_C1(25,2)/slope;
% strainOffset = exp12_C1(25,1) - yieldStrain;
% exp12_C1(1:25,1) = linspace(strainOffset,exp12_C1(25,1),25);
% for i = 1:25
%     exp12_C1(i,1) = (exp12_C1(i,2) - exp12_C1(1,2))/108000 + strainOffset;
% end
% exp12_C1(:,1) = exp12_C1(:,1) - exp12_C1(1,1);
% exp12_C1(:,1)=exp12_C1(:,1)-0.00178;
% save('exp12_C1.SS','exp12_C1','-ascii');
%%
% exp yield at about 241 at 0.0462 and 102 index
% mod reach ~242 at 0.0058 and 59 index
exp12_C2 = importdata('exp12_C2.SS');
exp12_C2x = exp12_C2(:,1)-0.0017;
exp12_C2y = exp12_C2(:,2);
mod12_C2 = importdata('mod12_C2.SS');
mod12_C2x = mod12_C2(:,1);
mod12_C2y = mod12_C2(:,2);

figure(2);
hold on;
plot(exp12_C2x,exp12_C2y,'--');
plot(mod12_C2x,mod12_C2y);
xlim([0,0.02])
%
% C2 12um exp data edit
% slope = mod12_C2(9,2)/mod12_C2(9,1);
% yieldStrain = exp12_C2(35,2)/slope;
% strainOffset = exp12_C2(35,1) - yieldStrain;
% for i = 1:35
%     exp12_C2(i,1) = (exp12_C2(i,2) - exp12_C2(1,2))/109350 + strainOffset;
% end
% exp12_C2(:,1) = exp12_C2(:,1) - exp12_C2(1,1);
% exp12_C2(:,1)=exp12_C2(:,1)-0.00178;


% save('exp12_C2.SS','exp12_C2','-ascii');
%%
% exp yield at about 382 at 0.0185 and 57 index
% mod reach ~382 at 0.0036 and 37 index
exp12_C3 = importdata('exp12_C3.SS');
exp12_C3x = exp12_C3(:,1);
exp12_C3y = exp12_C3(:,2);
mod12_C3 = importdata('mod12_C3.SS');
mod12_C3x = mod12_C3(:,1);
mod12_C3y = mod12_C3(:,2);

figure(3);
hold on;
plot(exp12_C3x,exp12_C3y,'--');
plot(mod12_C3x,mod12_C3y);
%
% C3 12um exp data edit
% slope = mod12_C3(3,2)/mod12_C3(3,1);
% yieldStrain = exp12_C3(34,2)/slope;
% strainOffset = exp12_C3(34,1) - yieldStrain;
% for i = 1:34
%     exp12_C3(i,1) = (exp12_C3(i,2) - exp12_C3(1,2))/126750 + strainOffset;
% end
% exp12_C3(:,1) = exp12_C3(:,1) - exp12_C3(1,1);
% save('exp12_C3.SS','exp12_C3','-ascii');

% exp yield at about 190 at 0.0247 and 38 index
% mod reach ~190 at 0.0036 and 37 index
exp20_C1 = importdata('exp20_C1.SS');
exp20_C1x = exp20_C1(:,1);
exp20_C1y = exp20_C1(:,2);
mod20_C1 = importdata('mod20_C1.SS');
mod20_C1x = mod20_C1(:,1);
mod20_C1y = mod20_C1(:,2);

figure(4);
hold on;
plot(exp20_C1x,exp20_C1y,'--');
plot(mod20_C1x,mod20_C1y);

% C1 20um exp data edit
% exp20_C1(:,1) = exp20_C1(:,1) - exp20_C1x(1);
% exp20_C1(:,2) = exp20_C1(:,2) - exp20_C1y(1);
% save('exp20_C1.SS','exp20_C1','-ascii');

% exp yield at about 251.2 at 0.0316 and 38 index
% mod reach ~250 at 0.0031 and 32 index
exp20_C2 = importdata('exp20_C2.SS');
exp20_C2x = exp20_C2(:,1);
exp20_C2y = exp20_C2(:,2);
mod20_C2 = importdata('mod20_C2.SS');
mod20_C2x = mod20_C2(:,1);
mod20_C2y = mod20_C2(:,2);

figure(5);
hold on;
plot(exp20_C2x,exp20_C2y,'--');
plot(mod20_C2x,mod20_C2y);

% C2 20um exp data edit
% exp20_C2(:,1) = exp20_C2(:,1) - exp20_C2x(1);
% exp20_C2(:,2) = exp20_C2(:,2) - exp20_C2y(1);
% save('exp20_C2.SS','exp20_C2','-ascii');

% exp yield at about 285 at 0.0233 and 45 index
% mod reach ~285 at 0.0033 and 34 index
exp20_C3 = importdata('exp20_C3.SS');
exp20_C3x = exp20_C3(:,1);
exp20_C3y = exp20_C3(:,2);
mod20_C3 = importdata('mod20_C3.SS');
mod20_C3x = mod20_C3(:,1);
mod20_C3y = mod20_C3(:,2);

figure(6);
hold on;
plot(exp20_C3x,exp20_C3y,'--');
plot(mod20_C3x,mod20_C3y);

% C3 20um exp data edit
% exp20_C3(:,1) = exp20_C3(:,1) - exp20_C3x(1);
% exp20_C3(:,2) = exp20_C3(:,2) - exp20_C3y(1);
% save('exp20_C3.SS','exp20_C3','-ascii');

% exp yield at about 229 at 0.0513 and 551 index
% mod reach ~280 at 0.0051 and 52 index
exp20_T1 = importdata('exp20_T1.SS');
exp20_T1x = exp20_T1(:,1);
exp20_T1y = exp20_T1(:,2);
mod20_T1 = importdata('mod20_T1.SS');
mod20_T1x = mod20_T1(:,1);
mod20_T1y = mod20_T1(:,2);

figure(7);
hold on;
plot(exp20_T1x,exp20_T1y,'--');
plot(mod20_T1x,mod20_T1y);

% exp yield at about 245.4 at 0.0225 and 226 index
% mod reach ~245 at 0.003 and 31 index
exp20_T2 = importdata('exp20_T2.SS');
exp20_T2x = exp20_T2(:,1);
exp20_T2y = exp20_T2(:,2);
mod20_T2 = importdata('mod20_T2.SS');
mod20_T2x = mod20_T2(:,1);
mod20_T2y = mod20_T2(:,2);

figure(8);
hold on;
plot(exp20_T2x,exp20_T2y,'--');
plot(mod20_T2x,mod20_T2y);

% exp yield at about 279 at 0.02686 and 298 index
% mod reach ~279 at 0.0031 and 31 index
exp20_T3 = importdata('exp20_T3.SS');
exp20_T3x = exp20_T3(:,1);
exp20_T3y = exp20_T3(:,2);
mod20_T3 = importdata('mod20_T3.SS');
mod20_T3x = mod20_T3(:,1);
mod20_T3y = mod20_T3(:,2);

figure(9);
hold on;
plot(exp20_T3x,exp20_T3y,'--');
plot(mod20_T3x,mod20_T3y);

% exp yield at about 208.7 at 0.03898 and 44 index
% mod reach ~209 at 0.0053 and 54 index
exp30_C2 = importdata('exp30_C2.SS');
exp30_C2x = exp30_C2(:,1);
exp30_C2y = exp30_C2(:,2);
mod30_C2 = importdata('mod30_C2.SS');
mod30_C2x = mod30_C2(:,1);
mod30_C2y = mod30_C2(:,2);

figure(10);
hold on;
plot(exp30_C2x,exp30_C2y,'--');
plot(mod30_C2x,mod30_C2y);

% C2 30um exp data edit
% slope = mod30_C2(5,2)/mod30_C2(5,1);
% yieldStrain = exp30_C2(1,2)/slope;
% strainOffset = exp30_C2(1,1) - yieldStrain;
% x = linspace(strainOffset,exp30_C2(1,1),20);
% y = linspace(0,exp30_C2(1,2),20);
% plot(x-x(1),y);
% exp30_C2 = [x',y';exp30_C2];
% exp30_C2(:,1) = exp30_C2(:,1) - exp30_C2(1,1);
% save('exp30_C2.SS','exp30_C2','-ascii');

% exp yield at about 285.5 at 0.001351 and 1 index
% mod reach ~286 at 0.0023 and 24 index
exp30_C3 = importdata('exp30_C3.SS');
exp30_C3x = exp30_C3(:,1);
exp30_C3y = exp30_C3(:,2);
mod30_C3 = importdata('mod30_C3.SS');
mod30_C3x = mod30_C3(:,1);
mod30_C3y = mod30_C3(:,2);

figure(11);
hold on;
plot(exp30_C3x,exp30_C3y,'--');
plot(mod30_C3x,mod30_C3y);

% C3 30um exp data edit
% slope = mod30_C3(2,2)/mod30_C3(2,1);
% yieldStrain = exp30_C3(1,2)/slope;
% strainOffset = exp30_C3(1,1) - yieldStrain;
% x = linspace(strainOffset,exp30_C3(1,1),20);
% y = linspace(0,exp30_C3(1,2),20);
% plot(x-x(1),y);
% exp30_C3 = [x',y';exp30_C3];
% exp30_C3(:,1) = exp30_C3(:,1) - exp30_C3(1,1);
% save('exp30_C3.SS','exp30_C3','-ascii');

%%
exp400_C1 = csvread('400-RD-1_TrueSSdata.csv');
exp400_C2= csvread('400-TD-1_TrueSSdata.csv');
exp400_C3= csvread('400-ND-1_TrueSSdata.csv');
exp625_C1 = csvread('625-RD-1_TrueSSdata.csv');
exp625_C2 = csvread('625-TD-1_TrueSSdata.csv');
exp625_C3 = csvread('625-ND-1_TrueSSdata.csv');

% 400C RD, yieldStrain = 0.0141, yieldStress = 181.91, approx. at smoothed
% data index 490
exp400_C1_smooth = smoothdata(exp400_C1(2:end,:),'movmean',10);
exp400_C1_smooth(:,1)=exp400_C1_smooth(:,1)-0.01579;

figure(12);
hold on;
plot(exp400_C1(:,1),exp400_C1(:,2));
plot([0;exp400_C1_smooth(:,1)],[0;exp400_C1_smooth(:,2)]);
linearx = 0:exp400_C1_smooth(1,1)/20:exp400_C1_smooth(1,1)-exp400_C1_smooth(1,1)/20;
lineary = 0:exp400_C1_smooth(1,2)/20:exp400_C1_smooth(1,2)-exp400_C1_smooth(1,2)/20;
exp400_C1_smooth = [linearx',lineary';exp400_C1_smooth];
% exp400_C1_smooth(:,1)=exp400_C1_smooth(:,1)-0.01579;

% figure;plot(exp400_C1_smooth(:,1),exp400_C1_smooth(:,2))
save('exp400_C1.SS','exp400_C1_smooth','-ascii');

A = importdata('exp400_C1.SS');
[x,y] = prepareCurveData(A(:,1),A(:,2));
[fit_o,~] = fit(x,y,'smoothingspline');
xin = linspace(0,x(end),1000);
plot(xin,fit_o(xin),'--');

% % Fit and find yield
% linearRegion = A(1:328,:); % row 328 pulled from examining the data
% % Fit to a line and get the slope
% linFit = polyfit(linearRegion(:,1),linearRegion(:,2),1);
% slope = linFit(1);
% intercept = linFit(2);
% % Make the 0.2% offset line
% linX = linspace(0,0.03,50);
% linY = @(x) slope*(x - 0.002);
% plot(linX,linY(linX));
% diff = @(x) fit_o(x) - linY(x);
% yieldStrain = fzero(diff,0);
% yieldStress = fit_o(yieldStrain);
%%
% 400C TD, yieldStrain = 0.0108, yieldStress = 263.8, approx. at smoothed
% data index 395
exp400_C2_smooth = smoothdata(exp400_C2(2:end,:),'movmean',10);
exp400_C2_smooth(:,1)=exp400_C2_smooth(:,1)-0.01192;
figure(13);
hold on;
plot(exp400_C2(:,1),exp400_C2(:,2));
plot([0;exp400_C2_smooth(:,1)],[0;exp400_C2_smooth(:,2)]);
linearx = 0:exp400_C2_smooth(1,1)/20:exp400_C2_smooth(1,1)-exp400_C2_smooth(1,1)/20;
lineary = 0:exp400_C2_smooth(1,2)/20:exp400_C2_smooth(1,2)-exp400_C2_smooth(1,2)/20;
exp400_C2_smooth = [linearx',lineary';exp400_C2_smooth];
save('exp400_C2.SS','exp400_C2_smooth','-ascii');

B = importdata('exp400_C2.SS');
[x,y] = prepareCurveData(B(:,1),B(:,2));
[fit_o,~] = fit(x,y,'smoothingspline');
xin = linspace(0,x(end),1000);
plot(xin,fit_o(xin),'--');

% % Fit and find yield
% % plot(B(:,1),B(:,2));
% linearRegion = B(1:152,:); % row 152 pulled from examining the data
% % Fit to a line and get the slope
% linFit = polyfit(linearRegion(:,1),linearRegion(:,2),1);
% slope = linFit(1);
% intercept = linFit(2);
% % Make the 0.2% offset line
% linX = linspace(0,0.03,50);
% linY = @(x) slope*(x - 0.002);
% plot(linX,linY(linX));
% diff = @(x) fit_o(x) - linY(x);
% yieldStrain = fzero(diff,0);
% yieldStress = fit_o(yieldStrain);

% 400C ND, yieldStrain = 0.0104, yieldStress = 230.99, approx. at smoothed
% data index 445
%%
exp400_C3_smooth = smoothdata(exp400_C3(2:end,:),'movmean',10);
exp400_C3_smooth(:,1)=exp400_C3_smooth(:,1)-0.01526;

figure(14);
hold on;
plot(exp400_C3(:,1),exp400_C3(:,2));
plot([0;exp400_C3_smooth(:,1)],[0;exp400_C3_smooth(:,2)]);
linearx = 0:exp400_C3_smooth(1,1)/20:exp400_C3_smooth(1,1)-exp400_C3_smooth(1,1)/20;
lineary = 0:exp400_C3_smooth(1,2)/20:exp400_C3_smooth(1,2)-exp400_C3_smooth(1,2)/20;
exp400_C3_smooth = [linearx',lineary';exp400_C3_smooth];
save('exp400_C3.SS','exp400_C3_smooth','-ascii');

C = importdata('exp400_C3.SS');
[x,y] = prepareCurveData(C(:,1),C(:,2));
[fit_o,~] = fit(x,y,'smoothingspline');
xin = linspace(0,x(end),1000);
plot(xin,fit_o(xin),'--');
%%
% % Fit and find yield
% % plot(C(:,1),C(:,2));
% linearRegion = C(1:271,:); % row 271 pulled from examining the data
% % Fit to a line and get the slope
% linFit = polyfit(linearRegion(:,1),linearRegion(:,2),1);
% slope = linFit(1);
% intercept = linFit(2);
% % Make the 0.2% offset line
% linX = linspace(0,0.03,50);
% linY = @(x) slope*(x - 0.002);
% plot(linX,linY(linX));
% diff = @(x) fit_o(x) - linY(x);
% yieldStrain = fzero(diff,0);
% yieldStress = fit_o(yieldStrain);

% 625C RD, yieldStrain = 0.0086, yieldStress = 173.04, approx. at smoothed
% data index 357
exp625_C1_smooth = smoothdata(exp625_C1(2:end,:),'movmean',10);
exp625_C1_smooth(:,1)=exp625_C1_smooth(:,1)-0.01;

figure(15);
hold on;
plot(exp625_C1(:,1),exp625_C1(:,2));
plot([0;exp625_C1_smooth(:,1)],[0;exp625_C1_smooth(:,2)]);
linearx = 0:exp625_C1_smooth(1,1)/20:exp625_C1_smooth(1,1)-exp625_C1_smooth(1,1)/20;
lineary = 0:exp625_C1_smooth(1,2)/20:exp625_C1_smooth(1,2)-exp625_C1_smooth(1,2)/20;
exp625_C1_smooth = [linearx',lineary';exp625_C1_smooth];
save('exp625_C1.SS','exp625_C1_smooth','-ascii');

D = importdata('exp625_C1.SS');
[x,y] = prepareCurveData(D(:,1),D(:,2));
[fit_o,~] = fit(x,y,'smoothingspline');
xin = linspace(0,x(end),1000);
plot(xin,fit_o(xin),'--');

% % Fit and find yield
% % plot(D(:,1),D(:,2));
% linearRegion = D(1:214,:); % row 214 pulled from examining the data
% % Fit to a line and get the slope
% linFit = polyfit(linearRegion(:,1),linearRegion(:,2),1);
% slope = linFit(1);
% intercept = linFit(2);
% % Make the 0.2% offset line
% linX = linspace(0,0.03,50);
% linY = @(x) slope*(x - 0.002);
% plot(linX,linY(linX));
% yieldStrain = fzero(@(x) fit_o(x) - linY(x),0.005);
% yieldStress = fit_o(yieldStrain);

% 625C TD, yieldStrain = 0.0103, yieldStress = 246.34, approx. at smoothed
% data index 268
exp625_C2_smooth = smoothdata(exp625_C2(2:end,:),'movmean',10);
exp625_C2_smooth(:,1)=exp625_C2_smooth(:,1)-0.015;
figure(16);
hold on;
plot(exp625_C2(:,1),exp625_C2(:,2));
plot([0;exp625_C2_smooth(:,1)],[0;exp625_C2_smooth(:,2)]);
linearx = 0:exp625_C2_smooth(1,1)/20:exp625_C2_smooth(1,1)-exp625_C2_smooth(1,1)/20;
lineary = 0:exp625_C2_smooth(1,2)/20:exp625_C2_smooth(1,2)-exp625_C2_smooth(1,2)/20;
exp625_C2_smooth = [linearx',lineary';exp625_C2_smooth];
save('exp625_C2.SS','exp625_C2_smooth','-ascii');

E = importdata('exp625_C2.SS');
[x,y] = prepareCurveData(E(:,1),E(:,2));
[fit_o,~] = fit(x,y,'smoothingspline');
xin = linspace(0,x(end),1000);
plot(xin,fit_o(xin),'--');

% % Fit and find yield
% % plot(E(:,1),E(:,2));
% linearRegion = E(1:134,:); % row 134 pulled from examining the data
% % Fit to a line and get the slope
% linFit = polyfit(linearRegion(:,1),linearRegion(:,2),1);
% slope = linFit(1);
% intercept = linFit(2);
% % Make the 0.2% offset line
% linX = linspace(0,0.03,50);
% linY = @(x) slope*(x - 0.002);
% plot(linX,linY(linX));
% yieldStrain = fzero(@(x) fit_o(x) - linY(x),0.005);
% yieldStress = fit_o(yieldStrain);

% 625C ND, yieldStrain = 0.0102, yieldStress = 218.49, approx. at smoothed
% data index 438
exp625_C3_smooth = smoothdata(exp625_C3(2:end,:),'movmean',15);
exp625_C3_smooth(:,1)=exp625_C3_smooth(:,1)-0.015;
figure(17);
hold on;
plot(exp625_C3(:,1),exp625_C3(:,2));
plot([0;exp625_C3_smooth(:,1)],[0;exp625_C3_smooth(:,2)]);
linearx = 0:exp625_C3_smooth(1,1)/20:exp625_C3_smooth(1,1)-exp625_C3_smooth(1,1)/20;
lineary = 0:exp625_C3_smooth(1,2)/20:exp625_C3_smooth(1,2)-exp625_C3_smooth(1,2)/20;
exp625_C3_smooth = [linearx',lineary';exp625_C3_smooth];
save('exp625_C3.SS','exp625_C3_smooth','-ascii');

F = importdata('exp625_C3.SS');
[x,y] = prepareCurveData(F(:,1),F(:,2));
[fit_o,~] = fit(x,y,'smoothingspline');
xin = linspace(0,x(end),1000);
plot(xin,fit_o(xin),'--');

% % Fit and find yield
% % plot(F(:,1),F(:,2));
% linearRegion = F(1:251,:); % row 251 pulled from examining the data
% % Fit to a line and get the slope
% linFit = polyfit(linearRegion(:,1),linearRegion(:,2),1);
% slope = linFit(1);
% intercept = linFit(2);
% % Make the 0.2% offset line
% linX = linspace(0,0.03,50);
% linY = @(x) slope*(x - 0.002);
% plot(linX,linY(linX));
% yieldStrain = fzero(@(x) fit_o(x) - linY(x),0.005);
% yieldStress = fit_o(yieldStrain);