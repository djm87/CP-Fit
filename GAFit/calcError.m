function [err,totN] = calcError(expX,expY,fitRange,simModx,simMody)

[simModelx, simModely] = prepareCurveData(simModx,simMody);
[simFitFcn, ~] = fit(simModelx,simModely,'smoothingspline');

if (fitRange(1) < fitRange(2))
    % limiting fitting region to a region with valid simulation and experiment
    % data
    startX = max([fitRange(1),expX(1),simModelx(1)]);
    iter = 1;
    while (expX(iter) < startX || simFitFcn(expX(iter))/expY(iter) < 0)
        iter = iter + 1;
    end
    startX = iter;
    
    endX = min([fitRange(2),expX(end),simModelx(end)]);
    iter = length(expX);
    while (expX(iter) > endX || simFitFcn(expX(iter))/expY(iter) < 0)
        iter = iter - 1;
    end
    endX = iter;
else
    startX = min([fitRange(1),expX(1),simModelx(1)]);
    iter = 1;
    while (expX(iter) > startX || simFitFcn(expX(iter))/expY(iter) < 0)
        iter = iter + 1;
    end
    startX = iter;

    endX = max([fitRange(2),expX(end),simModelx(end)]);
    iter = length(expX);
    while (expX(iter) < endX || simFitFcn(expX(iter))/expY(iter) < 0)
        iter = iter - 1;
    end
    endX = iter;
end

evalX = expX(startX:endX);

streak = 0;
for i = 1:100
    if (simFitFcn(evalX(i)) > 0)
        streak = streak + 1;
    else
        streak = 0;
    end
    if (streak >= 10)
        break;
    end
end
evalX = expX(startX+(i-10):endX);
evalY = expY(startX+(i-10):endX);

totN = length(evalX);

% error evaluation based on Tofallis, 2014
err = sqrt(sum((log(simFitFcn(evalX)./evalY)).^2)/totN);

end