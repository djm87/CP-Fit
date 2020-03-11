function info = readInfoFile(infoFname,cases)

%% Read from info file
file = fileread(infoFname);
fileText = regexp(file, '\r\n|\r|\n', 'split')';

% fitting strategy
info.fitStrat.ishift = str2num(fileText{3});
info.fitStrat.isave = str2double(fileText{5});
info.fitStrat.iplot = str2double(fileText{7});

% GA inputs
info.GAinp.gaType = str2double(fileText{10});
info.GAinp.iSingleRunGA = str2double(fileText{12});
if (fileText{14} == ' ')
    info.GAinp.A = [];
else
    info.GAinp.A = str2num(fileText{14});
end
if (fileText{16} == ' ')
    info.GAinp.b = [];
else
    info.GAinp.b = str2num(fileText{16});
end
if (fileText{18} == ' ')
    info.GAinp.Aeq = [];
else
    info.GAinp.Aeq = str2num(fileText{18});
end
if (fileText{20} == ' ')
    info.GAinp.beq = [];
else
    info.GAinp.beq = str2num(fileText{20});
end
if (fileText{22} == ' ')
    info.GAinp.nonlcon = [];
else
    info.GAinp.nonlcon = fileText{22};
end
info.GAinp.IntCon = str2num(fileText{24});

info.modelInfo.paramFileCount = str2num(fileText{27});

info.modelInfo.paramFileName = cell(info.modelInfo.paramFileCount,1);
for i = 1:info.modelInfo.paramFileCount
    info.modelInfo.paramFileName{i} = fileText{29+(i-1)};
end
skipline = info.modelInfo.paramFileCount-1;
if (skipline < 0)
    skipline = 0;
end

% system information
info.sysInfo.exeName = fileText{32+skipline};
info.sysInfo.slurmFlag = str2double(fileText{34+skipline});
info.sysInfo.parNode = str2double(fileText{36+skipline});
info.sysInfo.parTasks = str2double(fileText{38+skipline});
info.sysInfo.poolSize = str2double(fileText{40+skipline});
info.sysInfo.slurmMaxTime = fileText{42+skipline};
info.sysInfo.slurmJobName = fileText{44+skipline};
info.sysInfo.slurmPartition = fileText{46+skipline};
info.sysInfo.slurmMaxRam = str2double(fileText{48+skipline});
info.sysInfo.slurmNParallelJobs = str2double(fileText{50+skipline});

%% Check for additional fitting information and read the additional files
cyclicCaseIDs = cases{cases{1:end,{'IsCyclic'}} == true,{'CaseIdentifier'}};
info.cyclicFits = cell(length(cyclicCaseIDs),1);
for i = 1:length(cyclicCaseIDs)
    info.cyclicFits{i,1} = importdata(['FittingDataFiles/',cyclicCaseIDs{i},'.in']);
end

PFCaseIDs = cases{cases{1:end,{'IsPF'}} == true,{'CaseIdentifier'}};
info.PFFits = cell(length(PFCaseIDs),1);
for i = 1:length(cyclicCaseIDs)
    info.PFFits{i} = readPF(PFCaseIDs);
end
end