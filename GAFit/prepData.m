function expData = prepData(cases)

expData = [];

dataFiles = cases.FilePath;
X = cell(length(dataFiles),1);
smoothedY = cell(length(dataFiles),1);
expYFcn = cell(length(dataFiles),1);
parfor j = 1:length(dataFiles)
    importedExp = importVPSCout(dataFiles{j},0);
    X{j} = importedExp(:,1);
    expY = importedExp(:,2);
    smoothedY{j} = smooth(X{j},expY,0.005);
    
    [fitX, fitY] = prepareCurveData(X{j},smoothedY{j});
    [expYFcn{j}, ~] = fit(fitX,fitY,'smoothingspline');
end
expData.X = X;
expData.smoothedY = smoothedY;
expData.expYFcn = expYFcn;
end