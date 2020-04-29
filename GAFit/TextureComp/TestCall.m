%% Generate the random texture
tic
% Use MTEX to make a uniform ODF
CS = crystalSymmetry('622', [2.95 2.95 4.68], 'X||a', 'Y||b*', 'Z||c',...
    'mineral', 'Titanium (Alpha)', 'color', 'light blue');
SS=specimenSymmetry('-1');
odf=uniformODF(CS,SS);

% Export texture
fnameRand='random3000.TEX';
export_VPSC(odf,fnameRand,'points',3000);

%% Calc F for the random texture 
[euler,~]=loadVPSC(fnameRand);
euler=euler*degree;
FUniform = calcT(euler(:,1), euler(:,2), euler(:,3));

%% Calc F mean for target texture
fname='RD_EBSD_Initial.TEX'
[euler,weights]=loadVPSC(fname);
[FMeanTarget] = calcMeanF(euler*degree,weights);

%% Find the maximum difference between target texture and the random texture
maxdiff=max(sqrt(sum((FUniform'-FMeanTarget).^2)))

%% calc Error
fname='RD_EBSD_Initial.TEX'
[euler,weights]=loadVPSC(fname);
Fmean = calcMeanF(euler*degree,weights);
TDI = sqrt(sum((FMeanTarget-Fmean).^2))/maxdiff


fname='RD_EBSD_Initial_448.TEX';
[euler,weights]=loadVPSC(fname);
Fmean = calcMeanF(euler*degree,weights);
TDI = sqrt(sum((FMeanTarget-Fmean).^2))/maxdiff
toc
tic
fname='RD_EBSD_Initial_329.TEX';
[euler,weights]=loadVPSC(fname);
Fmean = calcMeanF(euler*degree,weights);
TDI = sqrt(sum((FMeanTarget-Fmean).^2))/maxdiff
toc
tic
fname='RD_EBSD_Initial_194.TEX';
[euler,weights]=loadVPSC(fname);
Fmean = calcMeanF(euler*degree,weights);
TDI = sqrt(sum((FMeanTarget-Fmean).^2))/maxdiff

fname='RD_EBSD_Initial_1.TEX';
[euler,weights]=loadVPSC(fname);
Fmean = calcMeanF(euler*degree,weights);
TDI = sqrt(sum((FMeanTarget-Fmean).^2))/maxdiff
toc