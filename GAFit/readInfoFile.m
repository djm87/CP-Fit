function info = readInfoFile(infoFname)

fileID = fopen(infoFname,'r');

fgetl(fileID); % System info title
fgetl(fileID); % System type
info.systemInfo.sysType = str2double(fgetl(fileID));
fgetl(fileID); % Program exe name
info.systemInfo.exeName = fgetl(fileID);
fgetl(fileID); % Section break

fgetl(fileID); % GA info title
fgetl(fileID); % Population number type
info.GAInfo.popNumType = str2double(fgetl(fileID));
fgetl(fileID); % Scaling factor
info.GAInfo.scaleFactor = str2double(fgetl(fileID));
fgetl(fileID); % Section break

fgetl(fileID); % Run parameter info title
fgetl(fileID); % Parameter table path
info.runInfo.paramTablePath = fgetl(fileID);
fgetl(fileID); % Number of headers in the parameter table
info.runInfo.paramNumTableHeaders = str2double(fgetl(fileID));
fgetl(fileID); % Table column headers
info.runInfo.paramTableHeaders = cell(1,info.runInfo.paramNumTableHeaders);
for i = 1:info.runInfo.paramNumTableHeaders
    info.runInfo.paramTableHeaders{i} = fgetl(fileID);
end
fgetl(fileID); % Recipe header
info.runInfo.paramRecipeHeader = fgetl(fileID);
fgetl(fileID); % Section break

fgetl(fileID); % Fitting strategy title
fgetl(fileID); % Shifting
info.fitInfo.ishift = str2double(fgetl(fileID));


fclose(fileID);

end