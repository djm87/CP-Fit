function yieldSearch(dataFileNames,systemParams)
    %% Shared variables breakdown
    % models = systemParams{1};
    rootPath = systemParams{2};
    % modelPaths = systemParams{3};
    cases = systemParams{4};
    % nPop = systemParams{5};
    % modelCaseDist = systemParams{6};
    expModNames = systemParams{7};
    % fitRange = systemParams{8};
    % iopt = systemParams{9};
    
    %% Instead of calculating the yield, the entire curve is being fitted.
    % The .TEX files will be transferred to respective optimization folders
    copied = 1;
    
    for i = 1:size(cases,1)
        curModeledCurves = cell(size(cases{i},1),1);
        for j = 1:size(cases{i},1)
            curModeledCurves{j} = importdata(dataFileNames{copied});
            copied = copied + 1;
        end
        save([rootPath,'/',expModNames{i}],'curModeledCurves');
    end
end