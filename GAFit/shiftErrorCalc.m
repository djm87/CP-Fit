function err = shiftErrorCalc(ishift,expX,expY,fitRange,simModx,simMody)

[simModelx, simModely] = prepareCurveData(simModx,simMody);
[simFitFcn, ~] = fit(simModelx,simModely,'smoothingspline');

if (sum(abs(ishift(1:2)))) % for some shifting amount of [-a,+b], -a < +b
    lb = ishift(1);
    ub = ishift(2);
    N = ishift(3);
    
    shiftErr = zeros(N,1);
    shiftX = linspace(lb,ub,N); 
    
    for shift = 1:N
        shiftAmount = shiftX(shift);
        tempx = expX + shiftAmount;
        
        % Fit original model/experiment to a function
        if length(tempx)>200
            [origDatax,origDatay] = prepareCurveData(downsample(tempx,200),...
                downsample(expY,200));
        else
            [origDatax,origDatay] = prepareCurveData(tempx,expY);                
        end
        [expFitFcn,~] = fit(origDatax,origDatay,'smoothingspline');
        
        adjMinFit = max([fitRange(1),tempx(1),simModelx(1)]);
        adjMaxFit = min([fitRange(2),tempx(end),simModelx(end)]);
        xIn = linspace(adjMinFit,adjMaxFit,100); % error calculation range

        funDiff = simFitFcn(xIn) - expFitFcn(xIn);
        shiftErr(shift) = sqrt(mean(funDiff.^2));
    end
    err = min(shiftErr);
    
else % if ishift is [0,0]
    % Fit original model/experiment to a function
    [origDatax,origDatay] = prepareCurveData(expX,expY);
    [expFitFcn,~] = fit(origDatax,origDatay,'smoothingspline');
    
    adjMinFit = max([fitRange(1),expX(1),simModelx(1)]);
    adjMaxFit = min([fitRange(2),expX(end),simModelx(end)]);
    xIn = linspace(adjMinFit,adjMaxFit,100);

    funDiff = weight*sqrt(sum((log(simFitFcn(xIn)./expFitFcn(xIn)))^2)/length(xIn));
    err = sqrt(mean(funDiff.^2));
end

end
