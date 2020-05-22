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
numconst = str2double(fileText{14});
if (numconst > 0)
    A_tot = [];
    b = zeros(numconst,1);
    for i = 1:numconst
        A = str2num(fileText{15+i});
        A_tot = [A_tot;A];
        b(i) = str2double(fileText{18+i});
    end
    info.GAinp.A = A_tot;
    info.GAinp.b = b;
else
    info.GAinp.A = [];
    info.GAinp.b = [];
end
offset = numconst;
if (fileText{20+offset} == ' ')
    info.GAinp.Aeq = [];
else
    info.GAinp.Aeq = str2num(fileText{20+offset});
end
if (fileText{22+offset} == ' ')
    info.GAinp.beq = [];
else
    info.GAinp.beq = str2num(fileText{22+offset});
end
if (fileText{24+offset} == ' ')
    info.GAinp.nonlcon = [];
else
    info.GAinp.nonlcon = fileText{24+offset};
end
info.GAinp.IntCon = str2num(fileText{26+offset});

info.modelInfo.paramFileCount = str2num(fileText{29+offset});

info.modelInfo.paramFileName = cell(info.modelInfo.paramFileCount,1);
for i = 1:info.modelInfo.paramFileCount
    info.modelInfo.paramFileName{i} = fileText{30+offset+i};
end
offset = offset + (info.modelInfo.paramFileCount-1);

% system information
info.sysInfo.exeName = fileText{34+offset};
info.sysInfo.slurmFlag = str2double(fileText{36+offset});
info.sysInfo.parNode = str2double(fileText{38+offset});
info.sysInfo.parTasks = str2double(fileText{40+offset});
info.sysInfo.poolSize = str2double(fileText{42+offset});
info.sysInfo.slurmMaxTime = fileText{44+offset};
info.sysInfo.slurmJobName = fileText{46+offset};
info.sysInfo.slurmPartition = fileText{48+offset};
info.sysInfo.slurmMaxRam = str2double(fileText{50+offset});
info.sysInfo.slurmNParallelJobs = str2double(fileText{52+offset});

%% Check for additional fitting information and read the additional files
cyclicCaseIDs = cases{cases{1:end,{'SS_Cyclic'}} == true,{'CaseIdentifier'}};
caseNo = cases{cases{1:end,{'SS_Cyclic'}} == true,{'CaseNo'}};
info.cyclicFits = cell(height(cases),1);
for i = 1:length(cyclicCaseIDs)
    info.cyclicFits{caseNo(i),1} = importVPSCout(['FittingDataFiles/',cyclicCaseIDs{i},'.in'],0);
end

end