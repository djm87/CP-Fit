function info = readInfoFile(infoFname)

file = fileread(infoFname);
fileText = regexp(file, '\r\n|\r|\n', 'split')';

% system information
info.sysInfo.sysType = str2double(fileText{3});
info.sysInfo.exeName = fileText{5};

% fitting strategy
info.fitStrat.ishift = str2num(fileText{8});
info.fitStrat.isave = str2double(fileText{10});
info.fitStrat.iplot = str2double(fileText{12});

% GA inputs
info.GAinp.gaType = str2double(fileText{15});
info.GAinp.iSingleRunGA = str2double(fileText{17});
if (fileText{19} == ' ')
    info.GAinp.A = [];
else
    info.GAinp.A = str2num(fileText{19});
end
if (fileText{21} == ' ')
    info.GAinp.b = [];
else
    info.GAinp.b = str2num(fileText{21});
end
if (fileText{23} == ' ')
    info.GAinp.Aeq = [];
else
    info.GAinp.Aeq = str2num(fileText{23});
end
if (fileText{25} == ' ')
    info.GAinp.beq = [];
else
    info.GAinp.beq = str2num(fileText{25});
end
if (fileText{27} == ' ')
    info.GAinp.nonlcon = [];
else
    info.GAinp.nonlcon = fileText{27};
end
info.GAinp.IntCon = str2num(fileText{29});

if (fileText{32} == ' ')
    info.modelInfo.paramFilePath = '/';
else
    info.modelInfo.paramFilePath = ['/',fileText{32},'/'];
end
info.modelInfo.paramFileName = fileText{34};

end