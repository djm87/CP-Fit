%Use: in matlab or octave navigate to the git repository and execute this script. 
%This installs the git repository to the startup paths for Matlab or Octave. 
%All scripts can then be executed from any directory in you system, allowing you to
%avoid having to copy the scripts around for different refinement directories.

P=genpath(pwd)
addpath(P)

%For root access use
%savepath

%For nonroot access specify a local directory and define the MATLABPATH to
%that directory. Matlab will then use the pathdef.m in the local directory.
savepath('/mnt/home/thrust2/djm87/pathdef.m')
