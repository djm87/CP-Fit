function [lc] = exportMAUDODFs(lc)
%exportODFs writes the ODFs in phase
    %extract the phase for the batch step
    phase=lc.phase{lc.BS};
    
    %Take care of the ODF directories
    udir=unique(lc.SamplePath(:,lc.BS));
    nudir=length(udir);
    for i=1:nudir
        odfsdir=fullfile(udir{i}, 'ODFs');
        odfWriteDir=fullfile(odfsdir,lc.OutputResult.Name{1,lc.BS}(1:end-4));
        %Make sure folders are present 
        if ~(exist(odfsdir,'dir')==7)
            mkdir(odfsdir);
        end
        if ~(exist(odfWriteDir,'dir')==7)
            mkdir(odfWriteDir);
        end
        lc = setODFPath(lc,udir{i},odfWriteDir);
    end                 
    for i=1:lc.ncases{lc.BS}
      if ~isempty(phase{i})
          for j=1:phase{i}.num %par file can have multiple phases
            if phase{i}.hasODF(j)==true
                lc.odfFileName{i,j,lc.BS}=fullfile(lc.odfWriteDir{i,lc.BS},...
                    [lc.CaseName{i,lc.BS} '_' phase{i}.name{j} '.MAUD']);
                writeMAUDODF(phase{i}.odf{j},lc.odfFileName{i,j,lc.BS})
            end
          end
      end
    end
end
function lc = setODFPath(lc,udir,odfWriteDir) 

 SP=lc.SamplePath(:,lc.BS);
 for j =1:lc.ncases{lc.BS}
    if strcmp(SP{j},udir)
         lc.odfWriteDir(j,lc.BS)={odfWriteDir};
    end
 end
end
function [] = writeMAUDODF(odf,odfFileName)
    fid=fopen(odfFileName,'w');
    for k=1:length(odf)
        fprintf(fid,'%8.4f %8.4f %8.4f %13.9f\n',...
            odf(k,1),odf(k,2),odf(k,3),odf(k,4));
    end
    fclose(fid);
end


