function GAoptions = readGAOptionsFile(GAFname)

file = fileread(GAFname);
fileText = regexp(file, '\r\n|\r|\n', 'split')';

totLines = size(fileText,1);
lineCount = 1;
while (lineCount <= totLines)
    if (fileText{lineCount}(1) == '*')
        lineCount = lineCount + 1;
    elseif (strcmp(fileText{lineCount},'gamultiobj') == 1 || strcmp(fileText{lineCount},'ga') == 1)
        % this part should only happen once at the start
        solver = fileText{lineCount};
        GAoptions = optimoptions(solver);
        lineCount = lineCount + 1;
    elseif (length(fileText{lineCount}) > 3 && strcmp(fileText{lineCount}(end-2:end),'Fcn') == 1)
        % next option is a function and may have function handles
        curOpt = fileText{lineCount};
        lineCount = lineCount + 1;
        fcnHdls = '{';
        while (fileText{lineCount}(1) == '@')
            fcnHdls = strcat(fcnHdls,fileText{lineCount},',');
            lineCount = lineCount + 1;
            if (fileText{lineCount}(1) == '>')
                if (isstrprop(fileText{lineCount}(2),'digit'))
                    fcnHdls = strcat(fcnHdls,fileText{lineCount}(2:end),',');
                    lineCount = lineCount + 1;
                else
                    fcnHdls = strcat(fcnHdls,'''',fileText{lineCount}(2:end),''',');
                    lineCount = lineCount + 1;
                end
            end
        end
        fcnHdls = strcat(fcnHdls(1:end-1),'}');
        GAoptions = optimoptions(GAoptions,curOpt,eval(fcnHdls));
    elseif (fileText{lineCount + 1}(1) == '+')
        curOpt = fileText{lineCount};
        optParam = fileText{lineCount + 1};
        if (strcmp(optParam(2:end),'true') == 1)
            optParam = true;
        elseif (isstrprop(optParam(2),'digit'))
            optParam = str2double(optParam(2:end));
        else
            optParam = optParam(2:end);
        end
        GAoptions = optimoptions(GAoptions,curOpt,optParam);
        lineCount = lineCount + 2;
    else
        curOpt = fileText{lineCount};
        GAoptions = optimoptions(GAoptions,curOpt);
        lineCount = lineCount + 1;
    end
end

end