function err = calcError(expX,expY,fitRange,simModx,simMody)

[simModelx, simModely] = prepareCurveData(simModx,simMody);
[simFitFcn, ~] = fit(simModelx,simModely,'smoothingspline');

% Fit original model/experiment to a function
[origDatax,origDatay] = prepareCurveData(expX,expY);
[expFitFcn,~] = fit(origDatax,origDatay,'smoothingspline');

adjMinFit = max([fitRange(1),expX(1),simModelx(1)]);
adjMaxFit = min([fitRange(2),expX(end),simModelx(end)]);
xIn = linspace(adjMinFit,adjMaxFit,1000);

% error evaluation based on Tofallis, 2014
err = sqrt(sum((log(simFitFcn(xIn)./expFitFcn(xIn)))^2)/length(xIn));

end