%Set environment
isOctave = exist('OCTAVE_VERSION', 'builtin') ~= 0;
if isOctave
    isMatlab=false;
else
    isMatlab=true;
end

%writeIns A script based generation of batch input file  
insName='WE43_refinement2.ins';

%Specify the sample directories
SampleDir={'//Split_Hopkinson_Bar/s1_2545_T5';
           '//Split_Hopkinson_Bar/s3_1185_T5';
           '//Split_Hopkinson_Bar/s4_1820';
           '//Split_Hopkinson_Bar/s5_1300_T5';
           '//Split_Hopkinson_Bar/s7_945';
           '//Split_Hopkinson_Bar/s9_770';
           '//Split_Hopkinson_Bar/s12_525_T5';};    
%Specify the folder names (must be an integer)
CaseNums={1; %Cases to run in each sample directory
          1;
          1;
          1;
          1;
          1;
          1};
assert(length(CaseNums)==length(SampleDir),'CaseNums needs to be same length as Sample directory')

%Specify the number of refinement steps
nSteps=1;
nStepsStart=3

%Generate the step names per sample, per refinement step
nSamples=length(CaseNums);
for i =1:nSamples
    for j=1:nSteps
        StepName{i,j}=['Step' int2str(nStepsStart+j)];
    end
end

%Specify the number of refinement itterations per sample, per
%refinement step
% NIter={5,5,3;
%        5,5,3}
NIter={};
NIterTmp=[4,4,4];
for i=1:nSamples
    for j=1:length(NIterTmp)
        NIter{i,j}=NIterTmp(j);
    end
end

%Specify the wizards per sample, per refinement step
% WizNum={1,13,14;
%         1,13,14};
WizNum={};
WizTmp =[6,13,14];
for i=1:nSamples
    for j=1:length(WizTmp)
        WizNum{i,j}=WizTmp(j);
    end
end

%Specify the input par per sample for the first refinement step
%Only the first case is needed.. others are generated below
% InputPar={'/initial.par'; 
%           '/initial.par'};
InputPar=cell(nSamples,nSteps);
InputParTmp='/Step3_Wiz14_Iters4.par';
for i=1:nSamples
    InputPar{i,1}=InputParTmp;
end


% Generate names per sample, per refinement step
% Uses OutputPar of previous refinement as InputPar for current
% refinement
for i=1:nSamples
    for j=1:nSteps
        OutputResult{i,j}=['/Batchresults_' StepName{i,j} '.txt']; %Only the starting case matters      
        OutputResultAutotrace{i,j}=['/Batchresults_Autotrace_' StepName{i,j} '.txt']; %These cases are specified in the par
        OutputPar{i,j}=['/' StepName{i,j} '_Wiz' int2str(WizNum{i,j}) '_Iters' int2str(NIter{i,j}) '.par'];
        if j>1
            InputPar{i,j}=OutputPar{i,j-1};
        end
    end
end

%Select from the available options (1==true, 0 ==false)
BatchDic = {"_riet_analysis_file"            ,  1;
           "_riet_analysis_iteration_number",  1;
           "_riet_analysis_wizard_index"    ,  1;
           "_riet_analysis_fileToSave"      ,  1;
           "_riet_meas_datafile_name"       ,  0;
           "_riet_append_simple_result_to"  ,  1;
           "_riet_append_result_to"         ,  1;
           "_riet_meas_datafile_replace"    ,  0;
           "_maud_background_add_automatic" ,  0;
           "_maud_output_plot_filename"     ,  0};
cnt=1;
for i=1:length(BatchDic)
    if BatchDic{i,2}
        BatchDicSelected(cnt)=i;
        cnt=cnt+1;
    end
end
BatchOptions=BatchDic(BatchDicSelected,1);

%open insName for writing
fid=fopen(insName,'w');

%loop over each refinement step and write to insName
for BS=1:nSteps

    %Write the CIF loop
    fprintf(fid,'loop_\n');
    for j=1:length(BatchOptions) 
     fprintf(fid,'%s\n',BatchOptions{j});
    end

    fprintf(fid,'\n');

    %For each sample write the cases to run
    for i=1:nSamples
        for j=1:length(CaseNums{i})
            id=CaseNums{i}(j);
            for k=1:length(BatchOptions)       
                switch BatchOptions{k}
                case "_riet_analysis_file"
                    BP=fullfile(SampleDir{i},int2str(id),InputPar{i,BS});
                    BP=checkBP(BP);
                    fprintf(fid,'''%s'' ',BP);
                case "_riet_analysis_wizard_index"
                    fprintf(fid,'%d ',WizNum{i,BS});
                case "_riet_analysis_iteration_number"
                    fprintf(fid,'%d ',NIter{i,BS});
                case "_riet_analysis_fileToSave"
                    BP=fullfile(SampleDir{i},int2str(id),OutputPar{i,BS});
                    BP=checkBP(BP);
                    fprintf(fid,'''%s'' ',BP);
                case "_riet_meas_datafile_name"
                  disp('Warning: using untested feature "_riet_meas_datafile_name"')  
                  %add when useful..
                  %BP=fullfile(SampleDir{i},int2str(id),measDatafileName{i,BS});
                  %fprintf(fid,'''%s'' ',BP);
                case "_riet_append_simple_result_to"
                   BP=fullfile(SampleDir{i},OutputResult{i,BS});
                   BP=checkBP(BP);
                   fprintf(fid,'''%s'' ',BP);
                case "_riet_append_result_to"
                   BP=fullfile(SampleDir{i},OutputResultAutotrace{i,BS});
                   BP=checkBP(BP);
                   fprintf(fid,'''%s'' ',BP);
                case "_riet_meas_datafile_replace"
                  disp('Warning: using untested feature "_riet_meas_datafile_replace"')  
                  %add when useful..
                  %fprintf(fid,'''%s'' ',MeasDatafileReplace{id,BS});
                case "_maud_background_add_automatic"
                  disp('Warning: using untested feature "_maud_background_add_automatic"')  
                  %add when useful
                  %fprintf(fid,'''%s'' ',AutoAddBK{id,BS});
                case "_maud_output_plot_filename"
                  disp('Warning: using untested feature "_maud_output_plot_filename"')  
                  %add when useful
                  %useless with multiple banks and rotations
                  %fprintf(fid,'''%s'' ',BatchPath{id,BS});
                end
            end %Switch
            
        %Go to next line
        fprintf(fid,'\n');
        
        end %Case loop
    end %Sample loop

    %Add space between refinement steps
    fprintf(fid,'\n');
end
    fclose(fid);

