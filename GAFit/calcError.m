function [err,totN] = calcError(expX,expY,fitRange,simModx,simMody)

[simModelx, simModely] = prepareCurveData(simModx,simMody);
[simFitFcn, ~] = fit(simModelx,simModely,'smoothingspline');

% limiting fitting region to a region with valid simulation and experiment
% data
startX = max(fitRange(1),expX(1),simModelx(1));
iter = 1;
while (expX(iter) < startX)
    iter = iter + 1;
end
startX = iter;

endX = min(fitRange(2),expX(end),simModelx(end));
iter = length(expX);
while (expX(iter) > endX)
    iter = iter - 1;
end
endX = iter;

evalX = expX(startX:endX);
evalY = expY(startX:endX);

totN = length(evalX);

% error evaluation based on Tofallis, 2014
err = sqrt(sum((log(simFitFcn(evalX)./evalY))^2)/totN);

end