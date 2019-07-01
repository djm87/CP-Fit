%% Initialize 
%Clear all variables
clear all 
tic
%Specify a file for batch 
lc.insName='WE43_refinement2.ins'
% lc.insName='WE43_refinement1.ins'

%% Read in. batch
lc=readIns(lc);

lc.options={ 'Max MAUD Instances',                           32; %used to set cpus or if slurm, tasks per node
             'Delete per cpu results',                        1;
             'Run MAUD Refinements',                          1;
             'Suppress cmd windows',                          0;%when true on PC Matlab does not wait for processes to end, but Octave works fine.
             'Use parallel for loops',                        0;%need parallel toolbox or octave parallel environement. Used in: custom refinements, phase extraction, odf operations
             'Using Slurm',                                   1;
             'Amount of mem per task (GB)',                   8; %slurm parameter Gb
             'Max ins',                                     500; %i.e. the number ins files 
             'Max slurm wall time (hr:min:sec)',      '01:00:00';
             'Run Title',                'WE43_Maud_Refinements';
             'Slurm Partion',                          'thrust2'}; %thrust2 or general

lc=setProgrammingEnvironment(lc);
         
for BS = 1:lc.numBatchStep  
    lc.BS = BS;
    lc = DistributeRun2CPU(lc);

    lc = WriteBatchInputs(lc);

    WriteCallParallel(lc); 

    runCases(lc);

    CleanRunDir;

    CleanupOuput(lc);
    
    lc = extractMAUDPhases(lc);
    
    lc = exportMAUDODFs(lc);   
    
    printMessage(sprintf('Done with Step %d\n\n',BS));
end
savelc(lc);
toc