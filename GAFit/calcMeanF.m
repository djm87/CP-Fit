function [meanF] = calcMeanF(euler,weights)
    %using the function made during "makeCalcTGSH.m"
    F = calcT(euler(:,1), euler(:,2), euler(:,3));
    %taking the weighted mean of the Fourier Coefficients
    meanF = (F')*weights./sum(weights);
end

