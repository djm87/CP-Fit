function [lc] = WriteBatchInputs(lc)
%WriteBatchInputs Writes the batch inputs for each cpu
    caseCnt=0;
    if lc.options{6,2}
        numParallelInstance=lc.ncases{lc.BS};
    else
        numParallelInstance=lc.options{1,2}
    end
    for i=1:numParallelInstance %Number of parallel instances
        insName=['Maud_Batch_input_' int2str(i) '.ins'];

        if ispc
            lc.batName{i}=['Maud_Batch_input_' int2str(i) '.bat'];
            fid=fopen(lc.batName{i},'w');
            fprintf(fid,'%s\n',...
             ['start /b /wait jre\bin\java -mx8192M -cp lib/Maud.jar;lib/ij.jar com.radiographema.MaudText -f %cd%\',insName]);
            fprintf(fid,'EXIT\n');
            fclose(fid);
        end
        
        fid=fopen(insName,'w');
        fprintf(fid,'loop_\n');
        for j=1:length(lc.BatchOptions(:,lc.BS)) 
            fprintf(fid,'%s\n',lc.BatchOptions{j,lc.BS});
        end
        
        fprintf(fid,'\n');
        
        for j=1:length(lc.caseID{i})
            id=lc.caseID{i}(j);

            for k=1:length(lc.BatchOptions(:,lc.BS))       
                switch lc.BatchOptions{k,lc.BS}
                    case '_riet_analysis_file'
                      fprintf(fid,'%s ',lc.InputPar.BatchPath{id,lc.BS});
                    case '_riet_analysis_wizard_index'
                      fprintf(fid,'%s ',lc.NIter{id,lc.BS});
                    case '_riet_analysis_iteration_number'
                      fprintf(fid,'%s ',lc.WizNum{id,lc.BS});
                    case '_riet_analysis_fileToSave'
                      fprintf(fid,'%s ',lc.OutputPar.BatchPath{id,lc.BS});
                    case '_riet_meas_datafile_name'
                      disp('Warning: using untested feature "_riet_meas_datafile_name"')  
                      %fprintf(fid,'%s ',lc.measDatafileName.BatchPath{id,lc.BS});
                    case '_riet_append_simple_result_to'
                      % handle potential write conflicts 
                      newstr=['_CPU' int2str(i) '.txt'];
                      lc.OutputResult.BatchPathCPU{id,lc.BS}=strrep(...
                          lc.OutputResult.BatchPath{id,lc.BS},'.txt',newstr);
                      lc.OutputResult.NameCPU{id,lc.BS}=strrep(...
                          lc.OutputResult.Name{id,lc.BS},'.txt',newstr);
                      fprintf(fid,'%s ',lc.OutputResult.BatchPathCPU{id,lc.BS});
                    case '_riet_append_result_to'
                      % handle potential write conflicts 
                      newstr=['_CPU' int2str(i) '.txt'];
                      lc.OutputResultAutotrace.BatchPathCPU{id,lc.BS}=strrep(...
                          lc.OutputResultAutotrace.BatchPath{id,lc.BS},'.txt',newstr);
                      lc.OutputResultAutotrace.NameCPU{id,lc.BS}=strrep(...
                          lc.OutputResultAutotrace.Name{id,lc.BS},'.txt',newstr);
                      fprintf(fid,'%s ',lc.OutputResultAutotrace.BatchPathCPU{id,lc.BS});
                    case '_riet_meas_datafile_replace'
                      disp('Warning: using untested feature "_riet_meas_datafile_replace"')  
                      %fprintf(fid,'%s ',lc.MeasDatafileReplace{id,lc.BS});
                    case '_maud_background_add_automatic'
                      disp('Warning: using untested feature "_maud_background_add_automatic"')  
                      %fprintf(fid,'%s ',lc.AutoAddBK{id,lc.BS});
                    case '_maud_output_plot_filename'
                      disp('Warning: using untested feature "_maud_output_plot_filename"')    
                      %fprintf(fid,'%s ',lc.OutputPlot.BatchPath{id,lc.BS});
                end
            end
            fprintf(fid,'\n');
                           
        end
        fclose(fid);
    end
end