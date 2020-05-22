function data = importVPSCout(path,nheaders)
% function reads VPSC output data and returns them as a m x n matrix

data = [];
fileID = fopen(path);
A = 0;
for i = 1:nheaders
    A = fgetl(fileID);
end
while (A ~= -1)
    A = fgetl(fileID);
    if (A ~= -1)
        data = [data;str2num(A)];
    end
end
fclose(fileID);

end