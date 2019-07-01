function setupParamFiles(caseID,populationSize,runSource,info)

paramPath = info.modelInfo.paramFilePath;
paramFname = info.modelInfo.paramFileName;

sourceFile = [runSource,'/',caseID,'/',paramFname];
for j = 1:populationSize
    curDir = ['RunningFolder/',caseID,'/',num2str(j)];
    destFile = [curDir,paramPath,paramFname];
    copyfile(sourceFile,destFile);
end

end