function [curSimData,errors] = errorEvalWrap(cases,nPop,caseIDs,runFoldName,curSimData,cyclicFits,errors,ishift)

FittingWeight = cases{1:end,'FittingWeight'};
fitRange = [cases.Start,cases.End];

outputFileName = cases.SimOut;
colX = cases.Column_x;
colY = cases.Column_y;
dataFiles = cases.FilePath;
iCyclic = cases.IsCyclic;
iPF = cases.IsPF;

for j = 1:numel(caseIDs)
    % Determine the case
    expData = importdata(dataFiles{j});
    expX = expData(:,1);
    expY = expData(:,2);
    curColX = colX(j);
    curColY = colY(j);
    curoutputFileName = outputFileName{j};
    curicyclic = iCyclic(j);
    curcyclicFits = cyclicFits{j};
    curFittingWeight = FittingWeight(j);
    curfitRange = fitRange(j,:);
    curiPF = iPF(j);
%     disp(j);
    parfor i = 1:nPop
%         disp(i);
        xdiff = 0; % shift data to start at (0,0) if not
        ydiff = 0;
        % save the data
        vpscData = importdata([runFoldName,'/',num2str(i),'/',num2str(j),'/',curoutputFileName]);
        vpscData = [vpscData.data(:,curColX),vpscData.data(:,curColY)];
        curSimData{i} = vpscData;
        
        % Calculate error ------------------------------------------------
        simModx = vpscData(:,1);
        simMody = vpscData(:,2);
        
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
                
                expY2_smoothed = smooth(expX2,expY2,0.1);
                
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

%                     plot(expX2,expY2_smoothed,'b','LineWidth',2);
%                     plot(simModx2,simMody2,'k','LineWidth',3);

                    if (ishift == 1)
                        [errorstmp(k),totSampletmp(k)] = calcErrorShift(expX2,expY2_smoothed,fitRange2,simModx2,simMody2);
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
            
        elseif (curiPF == 1)
            % Add pole figure error calculation here
%             if (info.fitStrat.ishift == 1)
%                 [errors(i,j),totSampletmp(k)] = FittingWeight(j)*calcErrorShift(expX,expY,fitRange(j,:),simModx,simMody);
%             else
%                 [errors(i,j),totSampletmp(k)] = FittingWeight(j)*calcError(expX,expY,fitRange(j,:),simModx,simMody);
%             end
        else
            if (ishift == 1)
                [errors(i,j),~] = calcErrorShift(expX,expY,curfitRange,simModx,simMody);
                if (~isreal(errors(i,j)))
                    errors(i,j) = 50;
                end
                errors(i,j) = curFittingWeight*errors(i,j);
            else
                [errors(i,j),~] = calcError(expX,expY,curfitRange,simModx,simMody);
                if (~isreal(errors(i,j)))
                    errors(i,j) = 50;
                end
                errors(i,j) = curFittingWeight*errors(i,j);
            end
        end
%         vpscData = importdata([runFoldName,'/',num2str(i),'/',num2str(j),'/ACT_PH1.OUT']);
%         vpscData = vpscData.data;
%         activitiesPH1{i} = vpscData;
%         vpscData = importdata([runFoldName,'/',num2str(i),'/',num2str(j),'/ACT_PH2.OUT']);
%         vpscData = vpscData.data;
%         activitiesPH2{i} = vpscData;
%         vpscData = importdata([runFoldName,'/',num2str(i),'/',num2str(j),'/ACT_PH3.OUT']);
%         vpscData = vpscData.data;
%         activitiesPH3{i} = vpscData;
    end
end

end