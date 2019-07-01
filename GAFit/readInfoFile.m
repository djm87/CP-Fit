function info = readInfoFile(infoFname)

file = fileread(infoFname);
fileText = regexp(file, '\r\n|\r|\n', 'split')';

% system information
info.sysInfo.sysType = str2double(fileText{3});
info.sysInfo.exeName = fileText{5};

% fitting strategy
info.fitStrat.ishift = str2double(fileText{8});

% GA inputs
info.GAinp.iSingleRunGA = str2double(fileText{11});
if (fileText{13} == ' ')
    info.GAinp.A = [];
else
    info.GAinp.A = str2num(fileText{13});
end
if (fileText{15} == ' ')
    info.GAinp.b = [];
else
    info.GAinp.b = str2num(fileText{15});
end
if (fileText{17} == ' ')
    info.GAinp.Aeq = [];
else
    info.GAinp.Aeq = str2num(fileText{17});
end
if (fileText{19} == ' ')
    info.GAinp.beq = [];
else
    info.GAinp.beq = str2num(fileText{19});
end
if (fileText{21} == ' ')
    info.GAinp.nonlcon = [];
else
    info.GAinp.nonlcon = fileText{21};
end
info.GAinp.IntCon = str2num(fileText{23});

end