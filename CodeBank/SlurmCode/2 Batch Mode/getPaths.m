function [lc]=getPaths(lc)
%getPaths Extracts paths needed for extracting the phase information 
    for i=1:lc.ncases
       [lc.path.samplePath{i},...
           lc.path.caseName{i},...
           lc.path.outputPar{i}]=GetSamplePath(lc.OutputPar{i});
       [lc.OutputResult(i).name,...
           lc.OutputResult(i).nameAfterCleanup]=...
           GetOutputResultDir(lc.OutputResult{i})
       [lc.path.OdfDir{i},...
            lc.path.OdfDirResult{i}]=GetODFDir(lc.path.samplePath{i},...
            lc.OutputResult{i}); 
    end
end

function [samplePath,caseName,outputParPath]=GetSamplePath(str2convert)
    str2convert=strrep(str2convert,'''','');
    str2convert=strrep(str2convert,'//','/');
    loc=strfind(str2convert,'/');
    caseName=str2convert(loc(end-1)+1:loc(end)-1);
    samplePath=str2convert(1:loc(end-1));
    if ispc 
       samplePath=strrep(samplePath,'/','\');
    end  
    samplePath=fullfile(pwd,samplePath);
    outputParPath=fullfile(pwd,str2convert);
end
function [name,nameAfterCleanup]=GetOutputResultDir(OutputResult)
    loc=strfind(OutputResult,'/');
    name=OutputResult(loc(end)+1:end-1)
    loc=strfind(output.name,'_');
    nameAfterCleanup=[output.name(1:loc(end)-1),'.txt']
    
end
function [OdfDir,OdfDirResult]=GetODFDir(samplePath,OutputResult)
    OdfDir=fullfile(samplePath,'ODF')
        loc=strfind(OutputResult,'/');

    
    OdfDirResult=fullfile(OdfDir,OutputResult(loc(end)+1:end-5));
end