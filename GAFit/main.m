function main()
%% Collecting information
% addpath(pwd);
% profile on;

inputFiles = fileread('Inp_GAFit_Ti.in');
inputFiles = regexp(inputFiles, '\r\n|\r|\n', 'split')';

caseFname = inputFiles{2};
fitParamFname = inputFiles{4};
recipeFname = inputFiles{6};
infoFname = inputFiles{8};
GAFname = inputFiles{10};
runSource = inputFiles{12};
runFoldName = inputFiles{14};

% Reading in the csv file containing information on the experimental data
% and where to find the simulation results
cases = readtable(caseFname);

% Reading in the csv file containing information on the parameters to be
% written and fitted for each generation. Also contains information on
% scaling, lb, ub for GA input, fitFlag, and more, see example file
fitParam = readtable(fitParamFname);

fitRecipe = readRecipeFile(recipeFname,size(fitParam,1));

% Reading in the text file containing various other information such as
% type of operating system and the name of the executable of the model to
% be run. ishift flag for considering a shifted curve in the errors.
% Also contains all the other GA inputs.
info = readInfoFile(infoFname,cases);

% Reading in the text file containing all the genetic algorithm options
% the user wishes to use. The function should be flexible enough to take on
% any options as long as various flags are used properly.
GAoptions = readGAOptionsFile(GAFname);

expData = prepData(cases);

% Set up pole figure calibration
PFFits = readPF(cases);

% Set up the run-time environment
% Creating a folder for each population and each run-case, load files are
% written at this time but parameter files will be written later when the
% population is generated
setupRunFolders(cases,GAoptions.PopulationSize,runSource,runFoldName);

%% Start optimization
optimizationStart(cases,fitParam,fitRecipe,info,GAoptions,runFoldName,PFFits,expData);

% profsave;

%% Delete the run folders to clean up space
% if (isunix)
%     system(['rm -r ',runFoldName]);
% else
%     system(['rmdir /S /Q ',runFoldName]);
% end

end