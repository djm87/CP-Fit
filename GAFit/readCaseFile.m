function cases = readCaseFile(caseFileName)

fileID = fopen(caseFileName,'r');

fgetl(fileID); % File header
fgetl(fileID); % Number of cases
cases.numCases = str2double(fgetl(fileID));
fgetl(fileID); % data file path explanation
fgetl(fileID); % data file path explanation cont

cases.caseList = cell(cases.numCases,1);
for i = 1:cases.numCases
    cases.caseList{i} = fgetl(fileID);
end

fclose(fileID);

end