function err = calcError(ishift,expX,expY,simFitFcn,fitRange)

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

        xIn = linspace(fitRange(1),fitRange(2),100); % error calculation range

        funDiff = simFitFcn(xIn) - expFitFcn(xIn);
        diffSq = funDiff.^2;
        diffSqSum = sum(diffSq);
        shiftErr(shift) = (diffSqSum);
    end
    err = min(shiftErr);
else % if ishift is [0,0]
    % Fit original model/experiment to a function
    [origDatax,origDatay] = prepareCurveData(expX,expY);
    [expFitFcn,~] = fit(origDatax,origDatay,'smoothingspline');
    
    xIn = linspace(fitRange(1),fitRange(2),100);

    funDiff = simFitFcn(xIn) - expFitFcn(xIn);
    diffSq = funDiff.^2;
    diffSqSum = sum(diffSq);
    err = (diffSqSum);
end

end
