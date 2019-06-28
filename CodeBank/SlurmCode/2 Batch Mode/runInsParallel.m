%% Initialize 
%Clear all variables
clear all 
tic
%Specify a file for batch 
lc.insName='Example.ins'

%% Read in. batch
lc=readIns(lc);

lc.options={'Max MAUD Instances',            4;
             'Delete per cpu results',       0;
             'Run MAUD Refinements',         0;
             'Suppress cmd windows',         0;%when true on PC Matlab does not wait for processes to end, but Octave works fine.
             'Use parallel for loops',       0};%need parallel toolbox or octave parallel environement. Used in: custom refinements, phase extraction, odf operations

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