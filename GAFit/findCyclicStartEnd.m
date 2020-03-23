function [simDataInc,simTmpXInds] = findCyclicStartEnd(simDataInc,simModx,simMody,expX2,expY2,iflip)
startInd = 0;
endInd = 0;
% end point of each cycle should be approximately the same, find the
% points in simulation data
for i = simDataInc:length(simModx)
    if (iflip*simModx(i) > iflip*expX2(1) && startInd == 0)
        if (simMody(i)/expY2(2) > 0)
            startInd = i;
        elseif (i < length(simModx) && iflip*simModx(i+1) < iflip*simModx(i))
            startInd = -1;
            break;
        end
    elseif (iflip*simModx(i) > iflip*expX2(end) && endInd == 0 || ...
            i == length(simModx) || (iflip*(simModx(i+1)) < iflip*(simModx(i))))
        endInd = i;
        break;
    end
end
simTmpXInds = [startInd,endInd];
simDataInc = i + 1;
end