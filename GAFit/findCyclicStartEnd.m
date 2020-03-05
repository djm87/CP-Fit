function [simDataInc,simTmpXInds] = findCyclicStartEnd(simDataInc,simModx,expX2,simTmpXInds,iflip)

if (iflip == 1)
    while (simDataInc ~= length(simModx) && (simModx(simDataInc+1) < simModx(simDataInc) || ...
            simModx(simDataInc+1) == simModx(simDataInc) ))
        simDataInc = simDataInc + 1;
    end
    while (simDataInc <= length(simModx))
        if (simModx(simDataInc) >= expX2(1) && simTmpXInds(1) == 0)
            simTmpXInds(1) = simDataInc;
        elseif ( ( (simModx(simDataInc) >= expX2(end)) || ...
                (simDataInc ~= 1 && (simModx(simDataInc) < simModx(simDataInc-1))) ) ...
                && simTmpXInds(2) == 0)
            if (simModx(simDataInc) < simModx(simDataInc-1))
                simTmpXInds(2) = simDataInc-1;
            else
                simTmpXInds(2) = simDataInc;
                simDataInc = simDataInc + 1;
            end
            break;
        elseif (simDataInc == length(simModx))
            simTmpXInds(2) = simDataInc;
        end
        simDataInc = simDataInc + 1;
    end
else
    while (simDataInc ~= 1 && (simModx(simDataInc+1) > simModx(simDataInc) || ...
                simModx(simDataInc+1) == simModx(simDataInc) ))
        simDataInc = simDataInc + 1;
    end
    while (simDataInc <= length(simModx))
        if (simModx(simDataInc) <= expX2(1) && simTmpXInds(1) == 0)
            simTmpXInds(1) = simDataInc;
        elseif ( ( (simModx(simDataInc) <= expX2(end)) || ...
                (simDataInc ~= 1 && (simModx(simDataInc) > simModx(simDataInc-1))) ) ...
                && simTmpXInds(2) == 0)
            if (simModx(simDataInc) < simModx(simDataInc-1))
                simTmpXInds(2) = simDataInc-1;
            else
                simTmpXInds(2) = simDataInc;
                simDataInc = simDataInc + 1;
            end
            break;
        elseif (simDataInc == length(simModx))
            simTmpXInds(2) = simDataInc;
        end
        simDataInc = simDataInc + 1;
    end
end

end