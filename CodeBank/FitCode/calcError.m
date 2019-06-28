function err = calcError(iopt,strainX,stressY,fitresult_mod,calcRange,varargin)
% Calculates the error depending on the iopt and various varargin would be
% required for various modes.
% The following inputs are always required in the listed order as they
% appear in almost all cases:
%       strainX and stressY - the original data to be compared to
%       fitresult_mod - current model's spline fit function handle
%       calcRange - range of data to compare
% 
% iopt = 1 calculates errors for various shifts of the specified range of
%       the curve and requires the following additional varargin inputs in
%       the order listed:
%       N - number of shifts
%       lb and ub - specifying the shifting bounds for the original data
% 
% iopt = 2 calculates error of the entire curve - model or not required
%       input is:
%       modDatax - the strain (x) data of the simulated model
%
% iopt = 3 calculates error of a portion of the curve (no shifting),
%       does not require any additional input.

switch (iopt)
    case 1
        if (length(varargin) ~= 3)
            error('Make sure the inputs are complete, see function header for details.');
        end
        N = varargin{1};
        lb = varargin{2};
        ub = varargin{3};
        
        % Calculates and returns errors of N shifts between given lower and upper
        % bounds lb, ub
        % the lowest error will be returned.
        shiftErr = zeros(N,1);
        shiftStrain = linspace(lb,ub,N);
        for shift = 1:N
            shiftAmount = shiftStrain(shift);
            tempx = strainX + shiftAmount;
            
            % Fit original model/experiment to a function
            if length(tempx)>200
            [origDatax,origDatay] = prepareCurveData(downsample(tempx,200),...
                downsample(stressY,200));
            else
            [origDatax,origDatay] = prepareCurveData(tempx,stressY);                
            end
            [fit_o,~] = fit(origDatax,origDatay,'smoothingspline');
            
            xIn = linspace(calcRange(1),calcRange(2),100); % error calculation range
            
            funDiff = fitresult_mod(xIn) - fit_o(xIn);
            diffSq = funDiff.^2;
            diffSqSum = sum(diffSq);
            shiftErr(shift) = (diffSqSum);
        end
        err = min(shiftErr);
    case 2
        if (length(varargin) ~= 1)
            error('Make sure the inputs are complete, see function header for details.');
        end
        modDatax = varargin{1};
        
        % Fit original model/experiment to a function
        [origDatax,origDatay] = prepareCurveData(strainX,stressY);
        [fit_o,~] = fit(origDatax,origDatay,'smoothingspline');
        
        % Take the maximum possible strain for fitting (in case
        % the model and data do not match)
        xIn = linspace(0,min(max(modDatax),max(origDatax)),100);
        
        funDiff = fitresult_mod(xIn) - fit_o(xIn);
        diffSq = funDiff.^2;
        diffSqSum = sum(diffSq);
        err = (diffSqSum);
        
    case 3
        % Fit original model/experiment to a function
        if length(strainX)>200
        [origDatax,origDatay] = prepareCurveData(downsample(strainX,200),...
            downsample(stressY,200));
        else
        [origDatax,origDatay] = prepareCurveData(strainX,stressY);                
        end
        [fit_o,~] = fit(origDatax,origDatay,'smoothingspline');
        
        % Take the maximum possible strain for fitting (in case
        % the model and data do not match)
        xIn = linspace(calcRange(1),calcRange(2),100);

        funDiff = fitresult_mod(xIn) - fit_o(xIn);
        diffSq = funDiff.^2;
        diffSqSum = sum(diffSq);
        err = (diffSqSum);
    otherwise
        error('iopt is not one of the valid options.');
end

end
