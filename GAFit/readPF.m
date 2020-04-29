function PFFits = readPF(cases)
% Calculates the F values for all the measured textures and stores the
% maximum difference as the difference between the reduced initial textures
% and the measured textures

PFFits = [];

% list of cases
PFIDs = cases{1:end,'DatasetIdentifier'};
PFRange = cases{1:end,'PFRange'};
caseIDs = cases{1:end,'CaseIdentifier'};

[uniquePFIDs,~] = unique(PFIDs,'stable');

[euler,~] = loadVPSC('FittingDataFiles/random3000.TEX',0);
euler=euler*pi/180;
FUniform = calcT(euler(:,1), euler(:,2), euler(:,3));

PFFits.maxDiff = zeros(length(uniquePFIDs),1);
% read and save initial textures and their F values
for i = 1:length(uniquePFIDs)
    [euler,weights] = loadVPSC(['FittingDataFiles/',uniquePFIDs{i},'_0.TEX'],0);
    FMeanTarget = calcMeanF(euler*pi/180,weights);
    % calculate max error as the error between measured initial textures and
    % one grain from a uniform texture that has the highest error
    maxdiff = max(sqrt(sum((FUniform'-FMeanTarget).^2)));
    
    PFFits.maxDiff(i) = maxdiff;
end

% expand the max diff to all the cases for ease of use later
maxDiff2 = zeros(length(caseIDs),1);
for i = 1:length(caseIDs)
    maxDiff2(i) = PFFits.maxDiff(cases{i,'Dataset'});
end
PFFits.maxDiff = maxDiff2;

PFFits.deformedF = cell(length(caseIDs),1);
PFFits.deformedMaxDiff = cell(length(caseIDs),1);
% read and save deformed textures and their F values
for i = 1:length(caseIDs)
    curRange = str2num(PFRange{i});
    
    PFFits.deformedF{i} = cell(length(curRange)-1,1);
    PFFits.deformedMaxDiff{i} = cell(length(curRange)-1,1);
    for j = 2:length(curRange)
        [euler,weights] = loadVPSC(['FittingDataFiles/',caseIDs{i},'_',num2str(curRange(j)*100),'.TEX'],0);
        FMeanTarget = calcMeanF(euler*pi/180,weights);
        PFFits.deformedF{i}{j-1} = FMeanTarget;
        PFFits.deformedMaxDiff{i}{j-1} = max(sqrt(sum((FUniform'-FMeanTarget).^2)));
    end
end

end