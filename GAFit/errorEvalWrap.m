function [curSimData,errors,shiftind] = errorEvalWrap(cases,nPop,runFoldName,...
    curSimData,errors,shiftind,PFFits,expData,cyclicFits,ishift)
% A wrapper that interacts with various functions to calculate errors for
% cyclic, regular stress-strain, and pole figure data

% Code currently supports the following objectives:
% 1 - Stress/strain error evaluation - option for cyclic
% 2 - Pole figure error evaluation
% 3 - Slope error evaluation
% 4 - Point error evaluation
totalObjectives = sum(cases.Objectives);

FittingWeight = cases{1:end,'FittingWeight'};
fitRange = [cases.SS_start,cases.SS_end];

caseIDs = cases{1:end,'CaseIdentifier'};
outputFileName = cases.SS_Fname;
colX = cases.SS_xCol;
colY = cases.SS_yCol;
dataFiles = cases.FilePath;
iCyclic = cases.SS_Cyclic;

for j = 1:numel(caseIDs)
    % Determine the case, regular cases are stress-strain curves, extended
    % cases are pole figures etc.
    expData.smoothedY = cell(length(dataFiles),1);
    expFcn = expData.expYFcn{j};
    expX = expData.X{j};
    expY = expData.smoothedY{j};
    curColX = colX(j);
    curColY = colY(j);
    curoutputFileName = outputFileName{j};
    curicyclic = iCyclic(j);
    curcyclicFits = cyclicFits{j};
    curFittingWeight = FittingWeight(j);
    curfitRange = fitRange(j,:);
    
    parfor i = 1:nPop

        xdiff = 0; % shift data to start at (0,0) if not
        ydiff = 0;
        % save the data
        vpscData = importVPSCout([runFoldName,'/',num2str(i),'/',num2str(j),'/',curoutputFileName],1);
        vpscData = [vpscData(:,curColX),vpscData(:,curColY)];
        curSimData{i,j} = vpscData;
%         vpscData = runData.lowestErrSimData{17,1}{j};
        % Calculate error ------------------------------------------------
        simModx = vpscData(:,1);
        simMody = vpscData(:,2);

        % cyclic error evaluation
        if (curicyclic == 1)
            % If cyclic, split the case into N cases. For each case,
            % calculate an error and take the average on the given indices
            errorstmp = zeros(size(curcyclicFits,1),1);
            totSampletmp = zeros(size(curcyclicFits,1),1);
            simDataInc = 1;
            iflip = 1;

            for k = 1:size(curcyclicFits,1)
                % Split experimental data
                expX2 = expX(curcyclicFits(k,1):curcyclicFits(k,2));
                expY2 = expY(curcyclicFits(k,1):curcyclicFits(k,2));

                expY2_smoothed = expY2;%smooth(expX2,expY2,0.005);

                if (k == 1)
                    xdiff = expX2(1);
                    ydiff = expY2_smoothed(1);
                    expX2 = expX2 - xdiff;
                    expY2_smoothed = expY2_smoothed - ydiff;
                else
                    expX2 = expX2 - xdiff;
                    expY2_smoothed = expY2_smoothed - ydiff;
                end

                % Find the section of the simulation data that matches
                if (iflip == 1)
                    [simDataInc,simTmpXInds] = findCyclicStartEnd(simDataInc,simModx,simMody,expX2,expY2,iflip);
                    iflip = -1;
                else
                    [simDataInc,simTmpXInds] = findCyclicStartEnd(simDataInc,simModx,simMody,expX2,expY2,iflip);
                    iflip = 1;
                end
                if (simTmpXInds(1) > 0 && simTmpXInds(2) > simTmpXInds(1) && simTmpXInds(2) - simTmpXInds(1) > 10)
                    simModx2 = simModx(simTmpXInds(1):simTmpXInds(2));
                    simMody2 = simMody(simTmpXInds(1):simTmpXInds(2));

                    fitRange2 = [expX2(1),expX2(end)];

                    if (ishift(1) == 1)
                        [errorstmp(k),totSampletmp(k),~] = calcErrorShift(expX2,expY2_smoothed,fitRange2,simModx2,simMody2,ishift);
                    else
                        [errorstmp(k),totSampletmp(k)] = calcError(expX2,expY2_smoothed,fitRange2,simModx2,simMody2);
                    end
                    if (~isreal(errorstmp(k)))
                        errorstmp(k) = 50;
                        totSampletmp(k) = 1;
                    end
                else
                    errorstmp(k) = 100;
                    totSampletmp(k) = 1;
                end
            end
            % errors are returned along with the number of samples (equal
            % to number of sampled data in experimental data. To calculate
            % average over all errors, the following mean is calculated:
            % error_i = sqrt( 1/N_i * sum( ln(sim/exp)^2 ) )
            % error = sum( error_i^2 * N_i / sum( N_i )
            errors(i,j) = curFittingWeight*sqrt( sum(errorstmp.^2.*totSampletmp) / sum(totSampletmp) );

        % non-cyclic stress-strain error evaluation
        else
            if (ishift(1) == 1)
                [errors(i,j),~,shiftind(i,j)] = calcErrorShift(expX,expFcn,curfitRange,simModx,simMody,ishift);
                if (~isreal(errors(i,j)))
                    errors(i,j) = 50;
                end
                errors(i,j) = curFittingWeight*errors(i,j);
            else
                [errors(i,j),~] = calcError(expX,expFcn,curfitRange,simModx,simMody);
                if (~isreal(errors(i,j)))
                    errors(i,j) = 50;
                end
                errors(i,j) = curFittingWeight*errors(i,j);
            end
        end % end of cyclic check
    end % end of parfor loop
end % end of case loop
    
iPF = cases.PFCount;
PFFileName = cases.PFFileName;
% evaluating texture error
if (any(iPF >= 1))
    iPFInd = find(iPF >= 1);
    texError = cell(nPop,numel(iPFInd));
    PFRange = cases.PFRange;
    totalPF = sum(cases{1:end,'PFCount'});
    for j = 1:numel(iPFInd)
        caseInd = iPFInd(j); % index of the corresponding case
        expF = PFFits.deformedF{caseInd};
        curMaxDiff = PFFits.maxDiff(caseInd);
        curPFFileName = PFFileName{caseInd};
        curPFRange = str2num(PFRange{caseInd});
        parfor i = 1:nPop
            % pole figure error evaluation
            texError{i,j} = zeros(1,length(curPFRange)-1);
            % fetch the simulation texture
            simTexture = [runFoldName,'/',num2str(i),'/',num2str(caseInd),'/',curPFFileName];
            % the first texture is also the texture at end of simulation
            if (curPFRange(2) == 100)
                [euler,weights] = loadVPSC(simTexture,100);
                FMeanTarget = calcMeanF(euler*pi/180,weights);
                
                diff = sqrt(sum((FMeanTarget-expF{1}).^2));
                texError{i,j}(1) = diff/curMaxDiff;
            else
                for k = 2:length(curPFRange)
                    [euler,weights] = loadVPSC(simTexture,curPFRange(k));
                    FMeanTarget = calcMeanF(euler*pi/180,weights);
                    diff = sqrt(sum((FMeanTarget-expF{k-1}).^2));
                    texError{i,j}(k-1) = diff/curMaxDiff;
                end
            end
        end
    end
    errors2 = zeros(nPop,size(errors,2)+totalPF);
    for i = 1:nPop
        errortmp = errors(i,:);
        for j = 1:numel(iPFInd)
            errortmp = [errortmp,texError{i,j}];
        end
        errors2(i,:) = errortmp;
    end
    errors = errors2;
end

end