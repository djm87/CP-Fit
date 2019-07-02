README/NOTES

When editing csv file, make sure to select all the rows to be deleted and use "delete rows" in Excel.
For some reason, if you simply use delete on the cells, MATLAB still sees those cells as part of the sheet.
It will attempt to read it in but it will replace all empty cells with '' or NaN, causing error.

Do not remove any lines starting with * and follow the format of the input values.
They are read in without any flexibility for simplicity of the code.
The only input that is flexible is the options for GA. In that case, be sure to follow the symbols.
If a new type of input is needed, feel free to contact me at zf1005@wildcats.unh.edu and I can update it.

All columns are required in the FittingParams.csv file and need to be named the same except for description column.

The program establishes the folders where each case and each population will be run once and for all at the beginning.
Make sure to place the respective files in the folders and make sure the folder names are specified correctly in the
Inp_WhatToFit.csv file.
The files that are common to all runs should be directly inside the source folder, currently EPSCSource.
The files that are unique to each case should be inside the respective subfolders inside the source folder.
The files that are unique to each population should be written directly into the run folders with WriteSxFile