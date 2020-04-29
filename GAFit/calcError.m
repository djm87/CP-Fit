function [err,totN] = calcError(expX,expFcn,fitRange,simModx,simMody)
% Calculates error for stress-strain curve or stress-strain curve sections

if (fitRange(1) < fitRange(2))
    % limiting fitting region to a region with valid simulation and experiment
    % data
    startX = max([fitRange(1),expX(1),simModx(1)]);
    iter = 1;
    while (simModx(iter) < startX || expFcn(expX(iter))/simMody(iter) < 0)
        iter = iter + 1;
    end
    startX = iter;
    
    endX = min([fitRange(2),expX(end),simModx(end)]);
    iter = length(simModx);
    while (simModx(iter) > endX || expFcn(expX(iter))/simMody(iter) < 0)
        iter = iter - 1;
    end
    endX = iter;
else
    startX = min([fitRange(1),expX(1),simModx(1)]);
    iter = 1;
    while (simModx(iter) > startX || expFcn(expX(iter))/simMody(iter) < 0)
        iter = iter + 1;
    end
    startX = iter;

    endX = max([fitRange(2),expX(end),simModx(end)]);
    iter = length(simModx);
    while (simModx(iter) < endX || expFcn(expX(iter))/simMody(iter) < 0)
        iter = iter - 1;
    end
    endX = iter;
end

evalX = simModx(startX:endX);
evalY = simMody(startX:endX);

totN = length(evalX);

% error evaluation based on Tofallis, 2014
err = sqrt(sum((log(evalY./expFcn(evalX))).^2)/totN);

end