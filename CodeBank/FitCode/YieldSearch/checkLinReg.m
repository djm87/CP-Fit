function [rsq,plStart] = checkLinReg(plStart,data)
    % plEnd depends on plStart, the range must cover at least 10% of the
    % total strain in the data, i.e. plEnd - plStart >= 10%(data(end,1) -
    % data(1,1))
    
    strainRange = 0.1*(data(end,1) - data(1,1));
    
    for i = plStart:length(data(:,1))
        if (data(i,1) > data(plStart,1)+strainRange)
            plEnd = i;
            break;
        end
    end
    
    xData = data(plStart:plEnd,1);
    yData = data(plStart:plEnd,2);
    p = polyfit(xData,yData,1);
    slope = p(1);

    yfit = polyval(p,xData);
    yresid = yData - yfit;
    ssresid = sum(yresid.^2);
    sstotal = (length(yData) - 1)*var(yData);
    rsq = 1 - ssresid/sstotal;

    plStart = plStart + 1;
end