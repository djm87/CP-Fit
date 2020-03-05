function [curSimData,errors] = errorEvalWrap(cases,nPop,caseIDs,runFoldName,curSimData,cyclicFits,errors,ishift)

FittingWeight = cases{1:end,'FittingWeight'};
fitRange = [cases.Start,cases.End];

outputFileName = cases.SimOut;
colX = cases.ColumnX;
colY = cases.ColumnY;
dataFiles = cases.FilePath;
iCyclic = cases.IsCyclic;
iPF = cases.IsPF;

for i = 1:nPop
    for j = 1:numel(caseIDs)
        % save the data
        vpscData = importdata([runFoldName,'/',num2str(i),'/',num2str(j),'/',outputFileName{j}]);
        vpscData = [vpscData.data(:,colX(j)),vpscData.data(:,colY(j))];
        curSimData{i} = vpscData;
        
        % Calculate error ------------------------------------------------
        simModx = vpscData(:,1);
        simMody = vpscData(:,2);
        
        % Determine the case
        expData = importdata(dataFiles{j});
        expX = expData(:,1);
        expY = expData(:,2);
        
        if (iCyclic(j) == 1)
            % If cyclic, split the case into N cases. For each case,
            % calculate an error and take the average on the given indices
            errorstmp = zeros(size(cyclicFits{j},1),1);
            totSampletmp = zeros(size(cyclicFits{j},1),1);
            simDataInc = 1;
            iflip = 1; % 1 looks for strain value greater, 2 looks for strain value lesser
            
            for k = 1:size(cyclicFits{j},1)
                % Split experimental data
                expX2 = expX(cyclicFits{j}(k,1):cyclicFits{j}(k,2));
                expY2 = expY(cyclicFits{j}(k,1):cyclicFits{j}(k,2));
                
                simTmpXInds = [0,0];
                % Find the section of the simulation data that matches
                if (iflip == 1)
                    [simDataInc,simTmpXInds] = findCyclicStartEnd(simDataInc,simModx,expX2,simTmpXInds,iflip);
                    iflip = 2;
                else
                    [simDataInc,simTmpXInds] = findCyclicStartEnd(simDataInc,simModx,expX2,simTmpXInds,iflip);
                    iflip = 1;
                end
                
                simModx2 = simModx(simTmpXInds(1):simTmpXInds(2));
                simMody2 = simMody(simTmpXInds(1):simTmpXInds(2));
                
                fitRange2 = [expX2(1),expX2(end)];
                
                if (ishift == 1)
                    [errorstmp(k),totSampletmp(k)] = calcErrorShift(expX2,expY2,fitRange2,simModx2,simMody2);
                else
                    [errorstmp(k),totSampletmp(k)] = calcError(expX2,expY2,fitRange2,simModx2,simMody2);
                end
            end
            % errors are returned along with the number of samples (equal
            % to number of sampled data in experimental data. To calculate
            % average over all errors, the following mean is calculated:
            % error_i = sqrt( 1/N_i * sum( ln(sim/exp)^2 ) )
            % error = sum( error_i^2 * N_i / sum( N_i )
            errors(i,j) = FittingWeight(j)*sqrt( sum(errorstmp.^2.*totSampletmp) / sum(totSampletmp) );
            
        elseif (iPF(j) == 1)
            % Add pole figure error calculation here
%             if (info.fitStrat.ishift == 1)
%                 [errors(i,j),totSampletmp(k)] = FittingWeight(j)*calcErrorShift(expX,expY,fitRange(j,:),simModx,simMody);
%             else
%                 [errors(i,j),totSampletmp(k)] = FittingWeight(j)*calcError(expX,expY,fitRange(j,:),simModx,simMody);
%             end
        else
            if (ishift == 1)
                [errors(i,j),~] = calcErrorShift(expX,expY,fitRange(j,:),simModx,simMody);
                errors(i,j) = FittingWeight(j)*errors(i,j);
            else
                [errors(i,j),~] = calcError(expX,expY,fitRange(j,:),simModx,simMody);
                errors(i,j) = FittingWeight(j)*errors(i,j);
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