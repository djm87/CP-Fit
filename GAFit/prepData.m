function expData = prepData(cases)

expData = [];

dataFiles = cases.FilePath;

expData.X = cell(length(dataFiles),1);
expData.smoothedY = cell(length(dataFiles),1);
expData.expYFcn = cell(length(dataFiles),1);
for j = 1:length(dataFiles)
    importedExp = importdata(dataFiles{j});
    expData.X{j} = importedExp(:,1);
    expY = importedExp(:,2);
    expData.smoothedY{j} = smooth(expData.X{j},expY,0.005);
    
    [fitX, fitY] = prepareCurveData(expData.X{j},expData.smoothedY{j});
    [expData.expYFcn{j}, ~] = fit(fitX,fitY,'smoothingspline');
end

end