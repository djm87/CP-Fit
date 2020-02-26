function writeDPSxFile(fname,inputs,caseID)
%Specify the case specific parameters using doubles, tripplets etc..

% fname is the full path that allows function to write the parameter file
% directly into each of the running folders
fileID = fopen(fname,'w');

% inputs has two columns - first column is all the parameters for this
% population for each parameter, second column is the scaling factor for the
% parameter.
par = inputs(:,1).*inputs(:,2);

% caseIDs provides information about the case that is currently writing.
% This can be used to distinguish the different grain size cases as long as
% it matches the given inputs.
% if (contains(caseID,'12'))
%     grainSize = 1.78*0.00001;
% elseif (contains(caseID,'20'))
%     grainSize = 2.22*0.00001;
% else
%     grainSize = 2.61*0.00001;
% end

fprintf(fileID,'Iron\n');
fprintf(fileID,'CUBIC             crysym\n');
fprintf(fileID,'   1.   1.   1.    90.   90.   90.   unit cell axes and angles\n');
fprintf(fileID,'Elastic stiffness for Fe at 300K [MPa] (Simmons and Huang)\n');
fprintf(fileID,' 206.37e3  135.0e3   135.0e3   000.0e3   000.0e3   000.0e3\n');
fprintf(fileID,' 135.0e3   206.37e3  135.0e3   000.0e3   000.0e3   000.0e3\n');
fprintf(fileID,' 135.0e3   135.0e3   206.37e3  000.0e3   000.0e3   000.0e3\n');
fprintf(fileID,' 000.0e3   000.0e3   000.0e3   117.0e3   000.0e3   000.0e3\n');
fprintf(fileID,' 000.0e3   000.0e3   000.0e3   000.0e3   117.0e3   000.0e3\n');
fprintf(fileID,' 000.0e3   000.0e3   000.0e3   000.0e3   000.0e3   117.0e3\n');
fprintf(fileID,'*Large elastic strain & pressure dependent Cij (kSM=0 or 1)\n');
fprintf(fileID,' 0 	\n');
fprintf(fileID,'*Thermal expansion coefficients of single crystal[K^(-1)]\n');
fprintf(fileID,'  5.7e-6   5.7e-6  10.3e-6   0.0e0   0.0e0   0.0e0\n');
fprintf(fileID,'SLIP AND TWINNING MODES\n');
fprintf(fileID,'3                               nmodesx\n');
fprintf(fileID,'2                               nmodes\n');
fprintf(fileID,'1  2                               mode(i)\n');
fprintf(fileID,'    {110}<111> SLIP\n');
fprintf(fileID,'  1   12   20    1   0                modex,nsmx,nrsx,isensex,itwx\n');
fprintf(fileID,'    0    1    1     1    1   -1       slip (n-b)\n');
fprintf(fileID,'    1    0    1     1    1   -1\n');
fprintf(fileID,'    1   -1    0     1    1   -1\n');
fprintf(fileID,'    0    1   -1     1   -1   -1\n');
fprintf(fileID,'    1    0    1     1   -1   -1\n');
fprintf(fileID,'    1    1    0     1   -1   -1\n');
fprintf(fileID,'    0    1    1     1   -1    1\n');
fprintf(fileID,'    1    0   -1     1   -1    1\n');
fprintf(fileID,'    1    1    0     1   -1    1\n');
fprintf(fileID,'    0    1   -1     1    1    1\n');
fprintf(fileID,'    1    0   -1     1    1    1\n');
fprintf(fileID,'    1   -1    0     1    1    1\n');
fprintf(fileID,'   {112}<111> SLIP\n');
fprintf(fileID,'  2   12   20    1   0                modex,nsmx,nrsx,isensex,itwx\n');
fprintf(fileID,'   -2    1   -1    -1   -1    1       slip (n-b)\n');
fprintf(fileID,'    1   -2   -1    -1   -1    1\n');
fprintf(fileID,'    1    1    2    -1   -1    1\n');
fprintf(fileID,'   -2   -1   -1    -1    1    1\n');
fprintf(fileID,'    1    2   -1    -1    1    1\n');
fprintf(fileID,'    1   -1    2    -1    1    1\n');
fprintf(fileID,'    2    1   -1     1   -1    1\n');
fprintf(fileID,'   -1   -2   -1     1   -1    1\n');
fprintf(fileID,'   -1    1    2     1   -1    1\n');
fprintf(fileID,'    2   -1   -1     1    1    1\n');
fprintf(fileID,'   -1    2   -1     1    1    1\n');
fprintf(fileID,'   -1   -1    2     1    1    1\n');
fprintf(fileID,'   {123}<111> SLIP\n');
fprintf(fileID,'  3   24   20    1   0                modex,nsmx,nrsx,isensex,itwx\n');
fprintf(fileID,'    1    2    3     1    1   -1       slip (n-b)\n');
fprintf(fileID,'   -1    3    2     1    1   -1\n');
fprintf(fileID,'    2    1    3     1    1   -1\n');
fprintf(fileID,'   -2    3    1     1    1   -1\n');
fprintf(fileID,'    3   -1    2     1    1   -1\n');
fprintf(fileID,'    3   -2    1     1    1   -1\n');
fprintf(fileID,'   -1    2   -3     1   -1   -1\n');
fprintf(fileID,'    1    3   -2     1   -1   -1\n');
fprintf(fileID,'    2   -1    3     1   -1   -1\n');
fprintf(fileID,'    2    3   -1     1   -1   -1\n');
fprintf(fileID,'    3    1    2     1   -1   -1\n');
fprintf(fileID,'    3    2    1     1   -1   -1\n');
fprintf(fileID,'    1   -2   -3     1   -1    1\n');
fprintf(fileID,'    1    3    2     1   -1    1\n');
fprintf(fileID,'    2   -1   -3     1   -1    1\n');
fprintf(fileID,'    2    3    1     1   -1    1\n');
fprintf(fileID,'    3    1   -2     1   -1    1\n');
fprintf(fileID,'    3    2   -1     1   -1    1\n');
fprintf(fileID,'    1    2   -3     1    1    1\n');
fprintf(fileID,'    1   -3    2     1    1    1\n');
fprintf(fileID,'    2    1   -3     1    1    1\n');
fprintf(fileID,'    2   -3    1     1    1    1\n');
fprintf(fileID,'   -3    1    2     1    1    1\n');
fprintf(fileID,'   -3    2    1     1    1    1\n');
fprintf(fileID,'DISLOCATION MODEL\n');
fprintf(fileID,' 0									!iDiag\n');
fprintf(fileID,'  0.9 0.0  25.0                       !INTERACTION CONSTANT, Q IN EQ. (3.14), grain size\n');
fprintf(fileID,'{110}<111> SLIP-----------------------------------------------\n');
fprintf(fileID,' 2.48000000e-10 9.00000000e-03\n');
fprintf(fileID,' %12.4e 7.00000000e+02\n',par(1));
fprintf(fileID,' 1.00000000e+07\n');
fprintf(fileID,' %12.4e 0.01\n',par(3));
fprintf(fileID,' %12.4e 0.00000000e+00 1.00000000e+00\n',par(2));
fprintf(fileID,' 0.   0.   0.                      !FOR HPFAC COEF FOR THIS SLIP MODE FOR GRAIN BOUNDARY, TWIN1 BOUNDARY, TWIN2 BOUNDARY\n');
fprintf(fileID,' 1.0 0. 30.                       !Q0,Q1,Q2 (K) FOR A IN EQ. (3.15) Q0+Q1*LOG(1+TEMP/Q2)\n');
fprintf(fileID,'{112}<111> SLIP-----------------------------------------------\n');
fprintf(fileID,' 2.48000000e-10 9.00000000e-03\n');
fprintf(fileID,' 2.50000000e+08 7.00000000e+02\n');
fprintf(fileID,' 1.00000000e+07\n');
fprintf(fileID,' 3.00000000e+11 0.01\n');
fprintf(fileID,' 4.50000000e+01 0.00000000e+00 1.00000000e+00\n');
fprintf(fileID,' 0.   0.   0.                      !FOR HPFAC COEF FOR THIS SLIP MODE FOR GRAIN BOUNDARY, TWIN1 BOUNDARY, TWIN2 BOUNDARY\n');
fprintf(fileID,' 1.0 0. 30.                       !Q0,Q1,Q2 (K) FOR A IN EQ. (3.15) Q0+Q1*LOG(1+TEMP/Q2)\n');
fprintf(fileID,' 1.0 0.5 0.1 0\n');
fprintf(fileID,' 0\n');
fprintf(fileID,' 1.0 0.95 0.95\n');
fprintf(fileID,' 1.0 0.95 0.95 \n');
fprintf(fileID,' 1.0 0.0 0.0\n');
fprintf(fileID,' 1.0 0.0 0.0 \n');
fprintf(fileID,'BACKSTRESS\n');
fprintf(fileID,'1.00000000e+02 1.00000000e+02       !tau_sat,ni\n');
fprintf(fileID,'3.00000000e-03 3.50000000e+00       !gam_b,fact\n');
fprintf(fileID,'1.00000000e+02 1.00000000e+02       !tau_sat,ni\n');
fprintf(fileID,'3.00000000e-03 3.50000000e+00       !gam_b,fact\n');
fprintf(fileID,'1\n');
fprintf(fileID,'1 0.0\n');
fclose(fileID);
end