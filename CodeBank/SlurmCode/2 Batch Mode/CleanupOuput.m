function [] = CleanupOuput(lc)
%CleanupOuput Combine result files into single file and delete per cpu results
%Note: script assumes that each sample has a single step to collapse for a
%given batch step (BS). To get around this constraint treat each case
%within the sample as a seperate sample with a unique batch result name.

    %Find the unique samples that were run
    udir=unique(lc.SamplePath(:,lc.BS));

    %Look through each sample directory
    for i =1:length(udir)
       
        %Read all cases into combinedResult combinedResultAutotrace
        combinedResult = readPartialResultFile(...
            lc.SamplePath(:,lc.BS),udir{i},...
            lc.OutputResult.NameCPU(:,lc.BS));
        combinedResultAutotrace = readPartialResultFile(...
            lc.SamplePath(:,lc.BS),udir{i},...
            lc.OutputResultAutotrace.NameCPU(:,lc.BS));
        
        %Sort the cases by gda number low to high
        ind = sortResultFile(combinedResult);
        
        %Only shuffled if ind is not empty
        %coded this way because gda not printed to autotrace file...
        combinedResult = reorderResultFile(combinedResult,ind);
        combinedResultAutotrace = reorderResultFile(combinedResultAutotrace,ind);

        %Write results
        ResultFilePath = getResultFileName(udir{i},...
            lc.OutputResult.Name{1,lc.BS});
        ResultAutotraceFilePath = getResultFileName(udir{i},...
            lc.OutputResultAutotrace.Name{1,lc.BS});
        writeCombinedResult(ResultFilePath,combinedResult);
        writeCombinedResult(ResultAutotraceFilePath,combinedResultAutotrace);

        removePerCoreFiles(ResultFilePath,lc.options{2,2});
        removePerCoreFiles(ResultAutotraceFilePath,lc.options{2,2});
    end
end
function [combinedResult] = readPartialResultFile(SP,udir,NameCPU)
% reads the per cpu result files into a combinedResult
    cntG=1;
    combinedResult=[];
    for j =1:length(SP)
        if strcmp(SP{j},udir)
            %File name to per cpu results
            resultFile=fullfile(udir,NameCPU{j});
            try
                %read file
                c=textread(resultFile,'%s','delimiter','\n');

                %Initialize the result file with the header for the first
                %per cpu result processed
                if cntG==1
                   combinedResult{cntG,1}=c{1}; 
                   cntG=cntG+1;
                end

                %Read in the case outputs
                cnt=2;
                while cnt<=length(c)
                    combinedResult{cntG,1}=c{cnt}; 
                    cnt=cnt+1;
                    cntG=cntG+1;
                end
            catch
                fprintf('File removed or did not run\n');
                fprintf('%s\n',resultFile);             
            end
        end
    end
end
function [ind] = sortResultFile(ResultFile)
%extracts the gda names from the result file and sorts them low to high
%returns the ordering
    loc=strfind(ResultFile(2:end),'.');
    gdaNum=zeros(length(loc),1);
    for j=1:length(loc)
        gdaNum(j)=str2num(ResultFile{j+1}(1:loc{j}(1)-1));
    end
    [~,ind]=sort(gdaNum); 
    ind=ind+1;
end
function [ResultFile] = reorderResultFile(ResultFile,ind)
%extracts reorders a result file according to ind
    if and(~isempty(ind),~isempty(ResultFile))
        tmp=ResultFile;
        for j=1:length(ind)
            ResultFile{j+1}=tmp{ind(j)};
        end 
    end
end
function [ResultFileName] = getResultFileName(udir,name)
% Creates the file name to write the combined results to   
    ResultFileName = fullfile(udir,name);
end
function [] = writeCombinedResult(filename,result)
%Writes the combined result of the CPUs
    if ~isempty(result)
        fid=fopen(filename,'w');
        for j=1:length(result)
            fprintf(fid,'%s\n',char(result(j)));
        end
        fclose(fid);   
    end
end
function [] = removePerCoreFiles(ResultPath,option)
%Removes all the partial result with CPU in name
    if option
        delete([ResultPath(1:end-4) '_CPU*']);
    end     
end
