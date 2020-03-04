function [err,totN] = calcErrorShift(expX,expY,fitRange,simModx,simMody)

lb = ishift(1);
ub = ishift(2);
N = ishift(3);

shiftErr = zeros(N,1);
totalN = zeros(N,1);
shiftX = linspace(lb,ub,N); 

for shift = 1:N
    shiftAmount = shiftX(shift);
    tempx = expX + shiftAmount;

    % Fit original model/experiment to a function
    if (length(tempx) > 200)
        [origDatax,origDatay] = prepareCurveData(downsample(tempx,200),...
            downsample(expY,200));
    else
        [origDatax,origDatay] = prepareCurveData(tempx,expY);                
    end
    
    [shiftErr(shift),totalN(shift)] = calcError(origDatax,origDatay,fitRange,simModx,simMody);
end
[err,ind] = min(shiftErr);
totN = totalN(ind);

end
