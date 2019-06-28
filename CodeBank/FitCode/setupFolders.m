function setupFolders(systemParams)
    %% Shared variables breakdown
    % models = systemParams{1};
    % rootPath = systemParams{2};
    modelPaths = systemParams{3};
    cases = systemParams{4};
    nPop = systemParams{5};
    modelCaseDist = systemParams{6};
    % curModel = systemParams{7};
    % weight = systemParams{8};
    % fitRange = systemParams{9};
    % iopt = systemParams{10};

    %% Function sets up all the synchronous folders according to the number
    % specified by maxsyn
    for i = 1:length(modelPaths)
        currentCases = modelCaseDist(i);
%         EPSCsource = [modelPaths{i},'/EPSC Source Files'];
        VPSCsource = [modelPaths{i},'/VPSC Source Files'];
        
        % make all the folders for all the models for all the cases
        for j = 1:currentCases
            for k = 1:nPop
                foldir = [modelPaths{i},'/',cases{i}(j,:),'_',int2str(k)];
                if (exist(foldir,'dir') ~= 7)
                    mkdir(modelPaths{i},[cases{i}(j,:),'_',int2str(k)]);
%                     copyfile(EPSCsource,foldir);
                    copyfile(VPSCsource,foldir);

                end
                
                filedir = [foldir,'/vpsc7.in'];
                if (exist(filedir,'file') ~= 2)
                    EPSC4FilePath = [modelPaths{i},'/vpsc7.in'];
                    copyfile(EPSC4FilePath,filedir);
                end
            end
        end
    end
end