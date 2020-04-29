function [err,totN,shiftind] = calcErrorShift(expX,expFcn,fitRange,simModx,simMody,ishift)

lb = ishift(2);
ub = ishift(3);
N = ishift(4);

shiftErr = zeros(N,1);
totalN = zeros(N,1);
shiftX = linspace(lb,ub,N); 

for shift = 1:N
    shiftAmount = shiftX(shift);
    tempx = simModx + shiftAmount;
    
    [shiftErr(shift),totalN(shift)] = calcError(expX,expFcn,fitRange,tempx,simMody);
end
[err,shiftind] = min(shiftErr);
totN = totalN(shiftind);

end
