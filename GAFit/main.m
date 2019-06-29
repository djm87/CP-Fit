function main()

fileID = fopen('InputFiles.in');
inputFiles = 
fclose(fileID);
caseFname = 'Inp_WhatToFit.csv';
infoFname = 'Inp_Info.in';
GAFname = 'Inp_GAOptions.in';

cases = readtable(caseFname);

info = readInfoFile(infoFname);

GAoptions = readGAOptionsFile(GAFname);

optimizationStart();

end