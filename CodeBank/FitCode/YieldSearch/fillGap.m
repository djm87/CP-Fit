function outArray = fillGap(cases,curModel,inArray)
    % This function converts each yield array into 8 elements with the
    % order of: C1/2/3, T1/2/3, PSC, SS for more uniform data processing
    space = ' ';
    checkCTPS = space(ones(1,8)); % total 8 possible loading scenarios
    % Check which loading conditions exist
    for i = 1:size(cases{curModel},1)
        curCase = cases{curModel}(i,1);
        if (curCase == 'C')
            if (cases{curModel}(i,2) == '1')
                checkCTPS(1) = '1';
            elseif (cases{curModel}(i,2) == '2')
                checkCTPS(2) = '1';
            elseif (cases{curModel}(i,2) == '3')
                checkCTPS(3) = '1';
            end
        elseif (curCase == 'T')
            if (cases{curModel}(i,2) == '1')
                checkCTPS(4) = '1';
            elseif (cases{curModel}(i,2) == '2')
                checkCTPS(5) = '1';
            elseif (cases{curModel}(i,2) == '3')
                checkCTPS(6) = '1';
            end
        elseif (curCase == 'P')
            checkCTPS(7) = '1';
        elseif (curCase == 'S')
            checkCTPS(8) = '1';
        end
    end
    
    % Since everything is in pre-determined order, there should be the same
    % number of yields as there are cases. Therefore, when the following
    % for loop finds the cases that exist (where checkCTPS is '1' not
    % 'NaN'), it will fill it up. Rest of the array would be zeros.
    outArray = zeros(1,8);
    curYield = 0;
    for i = 1:length(checkCTPS)
        if (~isnan(str2double(checkCTPS(i))))
            curYield = curYield + 1;
            outArray(i) = inArray(curYield);
        end
    end
end