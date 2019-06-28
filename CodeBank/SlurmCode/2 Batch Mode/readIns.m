function [lc] = readIns(lc)
%readIns Reads in .ins batch files into a struct
       
    c=textread(lc.insName,'%s','delimiter','\n');
    
    %Find the number of batch analyses
    loopStart=contains2(c,'loop_');
    loopLoc=find(loopStart);
    lc.numBatchStep=length(loopLoc);
    
    %read in each batch analysis
    for BS=1:lc.numBatchStep
        cnt=1;
        offset=loopLoc(BS);
        while ~isempty(c{cnt+offset})
            lc.BatchOptions{cnt,BS}=c{cnt+offset};
            cnt=cnt+1;
        end
        
        offset=cnt+offset;
        cnt=1;
        while ~isempty(c{cnt+offset})
            lc.cases2Run{cnt,BS}=c{cnt+offset};
            cnt=cnt+1;
        end
        lc.ncases{BS}=cnt-1;

        for i=1:lc.ncases{BS}
            curCase=strsplit(lc.cases2Run{i,BS});
            for j=1:length(lc.BatchOptions(:,BS))      
                switch lc.BatchOptions{j,BS}
                    case '_riet_analysis_file'
                      lc.InputPar.BatchPath{i,BS}=curCase{j};
                      lc.InputPar.FullPath{i,BS}=GetFullPath(curCase{j});
                      lc.InputPar.Name{i,BS}=GetName(lc.InputPar.FullPath{i,BS});
                    case '_riet_analysis_wizard_index'
                      lc.NIter{i,BS}=curCase{j};
                    case '_riet_analysis_iteration_number'
                      lc.WizNum{i,BS}=curCase{j}; 
                    case '_riet_analysis_fileToSave'
                      lc.OutputPar.BatchPath{i,BS}=curCase{j};
                      lc.OutputPar.FullPath{i,BS}=GetFullPath(curCase{j});
                      [lc.SamplePath{i,BS},lc.CaseName{i,BS}]=GetSamplePath(curCase{j});
                      lc.OutputPar.Name{i,BS}=GetName(lc.OutputPar.FullPath{i,BS});
                    case '_riet_meas_datafile_name'
                      disp('Warning: using untested feature "_riet_meas_datafile_name"')   
                      lc.measDatafileName.BatchPath{i,BS}=curCase{j};
                      lc.measDatafileName.FullPath{i,BS}=GetFullPath(curCase{j});
                      lc.measDatafileName.Name{i,BS}=GetName(lc.measDatafileName.FullPath{i,BS});
                    case '_riet_append_simple_result_to'
                      lc.OutputResult.BatchPath{i,BS}=curCase{j};
                      lc.OutputResult.FullPath{i,BS}=GetFullPath(curCase{j});
                      lc.OutputResult.Name{i,BS}=GetName(lc.OutputResult.FullPath{i,BS});
                    case '_riet_append_result_to'
                      lc.OutputResultAutotrace.BatchPath{i,BS}=curCase{j};
                      lc.OutputResultAutotrace.FullPath{i,BS}=GetFullPath(curCase{j});
                      lc.OutputResultAutotrace.Name{i,BS}=GetName(lc.OutputResultAutotrace.FullPath{i,BS});
                    case '_riet_meas_datafile_replace'
                      disp('Warning: using untested feature "_riet_meas_datafile_replace"')  
                      lc.MeasDatafileReplace{i,BS}=curCase{j};
                    case '_maud_background_add_automatic'
                      disp('Warning: using untested feature "_maud_background_add_automatic"')
                      lc.AutoAddBK{i,BS}=curCase{j};
                    case '_maud_output_plot_filename'
                      disp('Warning: using untested feature "_maud_output_plot_filename"')  
                      lc.OutputPlot.BatchPath{i,BS}=curCase{j};
                      lc.OutputPlot.FullPath{i,BS}=GetFullPath(curCase{j});
                      lc.OutputPlot.Name{i,BS}=GetName(lc.OutputPlot.FullPath{i,BS});
                end
            end
        end

    end
end
function [samplePath,caseName]=GetSamplePath(str2convert)
    str2convert=strrep(str2convert,'''','');
    if isunix
        FullPath=strrep(str2convert,'//','/');
        loc=strfind(FullPath,'/');
    elseif ispc
        FullPath=strrep(str2convert,'\\','\');
        loc=strfind(FullPath,'\');       
    end
    caseName=FullPath(loc(end-1)+1:loc(end)-1);
    samplePath=FullPath(1:loc(end-1));
    samplePath=fullfile(pwd,samplePath);
end
function [FullPath]=GetFullPath(str2convert)
    str2convert=strrep(str2convert,'''','');
    if isunix
        FullPath=strrep(str2convert,'//','/');
    elseif ispc 
        FullPath=strrep(str2convert,'\\','\');
    end
    FullPath=fullfile(pwd,FullPath);
end
function [Name]=GetName(str2convert)
    if ispc
        loc=strfind(str2convert,'\');
    elseif isunix
        loc=strfind(str2convert,'/');
    end
    Name=str2convert(loc(end)+1:end);

end