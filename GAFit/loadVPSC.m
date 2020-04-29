function [euler,weights] = loadVPSC(fname,strainref)
%loadVPSC loads the vpsc file and tries to delimiter options that are
%messing up the normal imputs

fileID = fopen(fname);
tline = fgetl(fileID);
while (ischar(tline))
    if (strcmp(tline(1:7),'TEXTURE') == 1)
        strain = str2double(tline(20:end));
        if (strain == strainref || strainref == 100)
            fgetl(fileID);
            fgetl(fileID);
            tline = fgetl(fileID);
            nGrain = str2double(tline(2:end));
            data = zeros(nGrain,4);
            for j = 1:nGrain
                tline = fgetl(fileID);
                data(j,:) = str2num(tline);
            end
            break;
        else
            tline = fgetl(fileID);
        end
    else
        tline = fgetl(fileID);
    end
end

fclose(fileID);

weights=data(:,4);
euler=data(:,1:3);

end