function WriteSxFile(fname,inputs,caseID)
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
if (strcmp(caseID(end-2:end),'12'))
    grainSize = 1;
elseif (strcmp(caseID(end-2:end),'20'))
    grainSize = 2;
else
    grainSize = 3;
end

fprintf(fileID,'*Material: Titanium                                                                                                               \n');
fprintf(fileID,'HEXAGONAL                                crysym                                                                                   \n');
fprintf(fileID,'1.  1.  1.587     90.  90.  120.         cdim(i),cang(i)                                                                          \n');
fprintf(fileID,'*Elastic stiffness of single crystal [MPa]                                                                                        \n');
fprintf(fileID,'  143.5e3    72.5e3    65.4e3     0.0     0.0     0.0                                                                             \n');
fprintf(fileID,'   72.5e3   143.5e3    65.4e3     0.0     0.0     0.0                                                                             \n');
fprintf(fileID,'   65.4e3    65.4e3   164.9e3     0.0     0.0     0.0                                                                             \n');
fprintf(fileID,'    0.0     0.0     0.0    32.1e3     0.0e3     0.0e3                                                                             \n');
fprintf(fileID,'    0.0     0.0     0.0     0.0e3    32.1e3     0.0e3                                                                             \n');
fprintf(fileID,'    0.0     0.0     0.0     0.0e3     0.0e3    35.5e3                                                                             \n');
fprintf(fileID,'*Thermal expansion coefficients of single crystal[K^(-1)] (currently for Zr)                                                      \n');
fprintf(fileID,'  5.7e-6   5.7e-6  10.3e-6   0.0e0   0.0e0   0.0e0                                                                                \n');
fprintf(fileID,'SLIP AND TWINNING MODES                                                                                                           \n');
fprintf(fileID,'8                               nmodesx                                                                                           \n');
fprintf(fileID,'5                               nmodes                                                                                            \n');
fprintf(fileID,'1  3 4 5 7                      mode(i)                                                                                           \n');
fprintf(fileID,'PRISMATIC                                                                                                                         \n');
fprintf(fileID,'  1    3   50   1                     modex,nsmx,nrsx,isensex                                                                     \n');
fprintf(fileID,'  0.0   0    0.   0.   0.   0.                twshx,isectw,thres1,thres2                                                          \n');
fprintf(fileID,'  45.5   42.  1290.   25.  75.0 0.    tau0x,tau1x,thet0,thet1, Hhp,Hgnd                                                           \n');
fprintf(fileID,'        1.0  1.0  10.0  2.0  0.0 0.0  0.0   hlatx(1,im,iph),im=1,nmodes                                                           \n');
fprintf(fileID,' 1  0 -1  0    -1  2 -1  0                                                                                                        \n');
fprintf(fileID,' 0 -1  1  0     2 -1 -1  0                                                                                                        \n');
fprintf(fileID,'-1  1  0  0    -1 -1  2  0                                                                                                        \n');
fprintf(fileID,'PYRAMIDAL<a>                                                                                                                      \n');
fprintf(fileID,'  2    6   50   1                     modex,nsmx,nrsx,isensex                                                                     \n');
fprintf(fileID,'  0.0   0    0.   0.   0.   0.                twshx,isectw,thres1,thres2                                                          \n');
fprintf(fileID,' 50.   50.  100.  100.   1.0  0.   tau0x,tau1x,thet0,thet1, Hhp,Hgnd                                                              \n');
fprintf(fileID,'        1.0  1.0  1.0  10.0 0.0 0.0  0.0   hlatx(1,im,iph),im=1,nmodes                                                            \n');
fprintf(fileID,' 1  0 -1  1    -1  2 -1  0                                                                                                        \n');
fprintf(fileID,' 0 -1  1  1     2 -1 -1  0                                                                                                        \n');
fprintf(fileID,'-1  1  0  1    -1 -1  2  0                                                                                                        \n');
fprintf(fileID,'-1  0  1  1    -1  2 -1  0                                                                                                        \n');
fprintf(fileID,' 0  1 -1  1     2 -1 -1  0                                                                                                        \n');
fprintf(fileID,' 1 -1  0  1     1  1 -2  0                                                                                                        \n');
fprintf(fileID,'PYRAMIDAL<c+a>                                                                                                                    \n');
fprintf(fileID,'  3   12   50   1                     modex,nsmx,nrsx,isensex                                                                     \n');
fprintf(fileID,'  0.0   0    0.   0.  0.   0.                twshx,isectw,thres1,thres2                                                           \n');
fprintf(fileID,'  495.  100.  1000.  5.  200.  0.      tau0x,tau1x,thet0,thet1, Hhp,Hgnd                                                          \n');
fprintf(fileID,'        1.0  1.0  2.0  2.0   0.0 0.0  0.0   hlatx(1,im,iph),im=1,nmodes                                                           \n');
fprintf(fileID,' 1  0 -1  1    -1 -1  2  3                                                                                                        \n');
fprintf(fileID,' 1  0 -1  1    -2  1  1  3                                                                                                        \n');
fprintf(fileID,' 0 -1  1  1     1  1 -2  3                                                                                                        \n');
fprintf(fileID,' 0 -1  1  1    -1  2 -1  3                                                                                                        \n');
fprintf(fileID,'-1  1  0  1     2 -1 -1  3                                                                                                        \n');
fprintf(fileID,'-1  1  0  1     1 -2  1  3                                                                                                        \n');
fprintf(fileID,'-1  0  1  1     2 -1 -1  3                                                                                                        \n');
fprintf(fileID,'-1  0  1  1     1  1 -2  3                                                                                                        \n');
fprintf(fileID,' 0  1 -1  1    -1 -1  2  3                                                                                                        \n');
fprintf(fileID,' 0  1 -1  1     1 -2  1  3                                                                                                        \n');
fprintf(fileID,' 1 -1  0  1    -2  1  1  3                                                                                                        \n');
fprintf(fileID,' 1 -1  0  1    -1  2 -1  3                                                                                                        \n');
fprintf(fileID,'BASAL                                                                                                                             \n');
fprintf(fileID,'  4    3   50   1                     modex,nsmx,nrsx,isensex                                                                     \n');
fprintf(fileID,'  0.0   0    0.   0.  0.   0.                 twshx,isectw,thres1,thres2                                                          \n');
fprintf(fileID,'  51.   52.  1011.  100.  1.  0.0    tau0x,tau1x,thet0,thet1, Hhp,Hgnd                                                            \n');
fprintf(fileID,'         1.0  1.0  1.0  10.0  0.0 0.0  0.0    hlatx(1,im,iph),im=1,nmodes                                                         \n');
fprintf(fileID,' 0  0  0  1     2 -1 -1  0                                                                                                        \n');
fprintf(fileID,' 0  0  0  1    -1  2 -1  0                                                                                                        \n');
fprintf(fileID,' 0  0  0  1    -1 -1  2  0                                                                                                        \n');
fprintf(fileID,'TENSILE TWINNING (TT1)                                                                                                            \n');
fprintf(fileID,'  5   6   50   0                     modex,nsmx,nrsx,isensex                                                                      \n');
fprintf(fileID,'  0.175   0    %5.1f  %5.1f  0.1  0.1         twshx,isectw,thres1,thres2     Note thres1 is the inverse of the # of lamella          \n',par(36),par(37));
fprintf(fileID,'  165.  20.  100.  0. 135.  0. 0.0 0.0    tau0x,tau1x,thet0,thet1, Hhp,Hgnd                                                       \n');
fprintf(fileID,'         1.0  1.0   10.0  15.0 1.0 1.0 0.0 hlatx(1,im,iph),im=1,nmodes                                                            \n');
fprintf(fileID,' 1  0 -1  2    -1  0  1  1                                                                                                        \n');
fprintf(fileID,' 0  1 -1  2     0 -1  1  1                                                                                                        \n');
fprintf(fileID,'-1  1  0  2     1 -1  0  1                                                                                                        \n');
fprintf(fileID,'-1  0  1  2     1  0 -1  1                                                                                                        \n');
fprintf(fileID,' 0 -1  1  2     0  1 -1  1                                                                                                        \n');
fprintf(fileID,' 1 -1  0  2    -1  1  0  1                                                                                                        \n');
fprintf(fileID,'TENSILE TWINNING (TT2)                                                                                                            \n');
fprintf(fileID,'  6    6   50   0                     modex,nsmx,nrsx,isensex                                                                     \n');
fprintf(fileID,'  0.6277   1    0.2  0.5   0.2  0.5             twshx,isectw,thres1,thres2   Note thres1 is the inverse of the # of lamella       \n'); %not sure if # of lamella is appropriate
fprintf(fileID,'  340.   340.  1.  0.50  100.  0.  0.0 0.0    tau0x,tau1x,thet0,thet1, Hhp,Hgnd                                                   \n');
fprintf(fileID,'        1.0  1.0   10.0   5.0 1.0 1.0  0.0 hlatx(1,im,iph),im=1,nmodes                                                            \n');
fprintf(fileID,' 1  1 -2  1    1  1 -2 -6                                                                                                         \n');
fprintf(fileID,'-1  2 -1  1   -1  2 -1 -6                                                                                                         \n');
fprintf(fileID,'-2  1  1  1   -2  1  1 -6                                                                                                         \n');
fprintf(fileID,'-1 -1  2  1   -1 -1  2 -6                                                                                                         \n');
fprintf(fileID,' 1 -2  1  1    1 -2  1 -6                                                                                                         \n');
fprintf(fileID,' 2 -1 -1  1    2 -1 -1 -6                                                                                                         \n');
fprintf(fileID,'COMPRESSIVE TWINNING (CT1)                                                                                                        \n');
fprintf(fileID,'  7    6   50   0                     modex,nsmx,nrsx,isensex                                                                     \n');
fprintf(fileID,'  0.218   0    %5.1f  %5.1f   0.1  0.4             twshx,isectw,thres1,thres2   Note thres1 is the inverse of the # of lamella        \n',par(40),par(41)); %not sure if # of lamella is appropriate
fprintf(fileID,'  340.   340.  1.  0.50  100.  0.  0.0 0.0   tau0x,tau1x,thet0,thet1, Hhp,Hgnd                                                    \n');
fprintf(fileID,'        1.0  1.0   10.0   5.0 1.0  1.0  0.0    hlatx(1,im,iph),im=1,nmodes                                                        \n');
fprintf(fileID,' 2 -1 -1  2     2 -1 -1 -3                                                                                                        \n');
fprintf(fileID,' 1  1 -2  2     1  1 -2 -3                                                                                                        \n');
fprintf(fileID,'-1  2 -1  2    -1  2 -1 -3                                                                                                        \n');
fprintf(fileID,'-2  1  1  2    -2  1  1 -3                                                                                                        \n');
fprintf(fileID,'-1 -1  2  2    -1 -1  2 -3                                                                                                        \n');
fprintf(fileID,' 1 -2  1  2     1 -2  1 -3                                                                                                        \n');
fprintf(fileID,'COMPRESSIVE TWINNING (CT2)                                                                                                        \n');
fprintf(fileID,'  8    6   50   0                     modex,nsmx,nrsx,isensex                                                                     \n');
fprintf(fileID,'  0.1043   1    0.2  0.9   0.2  0.5             twshx,isectw,thres1,thres2   Note thres1 is the inverse of the # of lamella       \n'); %not sure if # of lamella is appropriate
fprintf(fileID,'  340.   340.  1.  0.50  100.  0.  0.0 0.0   tau0x,tau1x,thet0,thet1, Hhp,Hgnd                                                    \n');
fprintf(fileID,'        1.0  1.0   10.0   5.0  1.0  1.0  0.0   hlatx(1,im,iph),im=1,nmodes                                                        \n');
fprintf(fileID,' 1  0 -1  1     1  0 -1 -2                                                                                                        \n');
fprintf(fileID,' 0  1 -1  1     0  1 -1 -2                                                                                                        \n');
fprintf(fileID,'-1  1  0  1    -1  1  0 -2                                                                                                        \n');
fprintf(fileID,'-1  0  1  1    -1  0  1 -2                                                                                                        \n');
fprintf(fileID,' 0 -1  1  1     0 -1  1 -2                                                                                                        \n');
fprintf(fileID,' 1 -1  0  1     1 -1  0 -2                                                                                                        \n');
fprintf(fileID,'-------------------                                                                                                               \n');
fprintf(fileID,'DISLOCATION MODEL                                                                                                                 \n');
fprintf(fileID,'  %9.7g 1.0 1                       !GRSZE, SCLGAMD0, iddsys                                                                       \n',par(grainSize)); %mean grain size used in HP equations for twins grain size set by aspect ratios in .in
fprintf(fileID,'  0.086 1.380622e-23 1.e-03        !k_deb, BOLTZMAN (J/K),  REF STRAIN RATE (1/s)   k_deb hard coded in EPSC see 3.20             \n');
fprintf(fileID,'PRISMATIC                                                                                                                         \n');
oset=6;
fprintf(fileID,' 2.95e-10 0.53779 %7.5g            !BURG (m), KOVB3 (J/Km^3), NORM ACTENER g IN EQ. (3.12)                                        \n',par(oset+3));
fprintf(fileID,' %7.5g %7.5g                  !KGENER-K1 IN EQ. (3.8) (1/m), DRAG STRESS-D IN EQ. (3.12) (Pa)                                \n',par(oset+2),par(oset+4));
fprintf(fileID,' 1.E+7                             ! EDOT_O IN EQ. (3.12)                                                                         \n');
fprintf(fileID,' %9.7g 1.E-01                    !INITIAL RHO_S (1/m^2), INITIAL RHO_DEB FOR EACH SLIP MODE (1/m^2)                             \n',par(oset+7));
fprintf(fileID,' 0 %5.1f   0.   1.0   0.  0. 0.  0.0  0.0	  !ISL (law for slip),A,B,C,D,E,... SLIP CONSTANTS (set for no temperature dependence)\n',par(oset)); %This is not a function of temperature. Our data is not a function of temperature so probably ok.
fprintf(fileID,' 2   %9.7g  %9.7g                 !TLATENT HARDENING BY THIS SLIP MODE ON TWIN1, TWIN2:  C IN EQ. (3.28                          \n',par(oset+5),par(oset+6));
fprintf(fileID,' 0.  0.                            !TLATENT1 HARDENING BY THIS SLIP MODE ON TWIN1, TWIN2                                          \n');
fprintf(fileID,' %5.1f  %5.1f %5.1f                     !FOR HPFAC COEF FOR THIS SLIP MODE FOR GRAIN BOUNDARY, TWIN1 BOUNDARY, TWIN2 BOUNDARY          \n',par(oset+1),par(oset+1),par(oset+1)); %Going to try hitting the same HP for twin and grain.
fprintf(fileID,' 2 8.0 0. 1.                       !IQ (LAW FOR Q), Q0,Q1,Q2 FOR LAW FOR Q, Q2 just needs to be non-zero                          \n');
fprintf(fileID,' 2 0.0 0.0 0.0 44.0e+03            !ISHMOD (law for shear modulus),SHEAR MODULUS CONSTANTS                                        \n');oset=24;
fprintf(fileID,'PYRAMIDAL<c+a>                                                                                                                    \n');
fprintf(fileID,' 5.54e-10  8.1198-02 %7.5g         !BURG (m), KOVB3 (J/Km^3), NORM ACTENER g IN EQ. (3.12)                                         \n',par(oset+3));
fprintf(fileID,' %7.5g  %7.5g                 !KGENER-K1 IN EQ. (3.8) (1/m), DRAG STRESS-D IN EQ. (3.12) (Pa)                                \n',par(oset+2),par(oset+4));
fprintf(fileID,' 1.E+07                            !EDOT_O IN EQ. (3.12)                                                                          \n');
fprintf(fileID,' %9.7g 1.E-01                    !INITIAL RHO_S (1/m^2), INITIAL RHO_DEB FOR EACH SLIP MODE (1/m^2)                            \n',par(oset+7));
fprintf(fileID,' 0 %5.1f 0.  1.  0.  0. 	0. 0.0 0.0   !ISL (law for slip),A,B,C,D,E,... SLIP CONSTANTS                                             \n',par(oset));
fprintf(fileID,' 2 %9.7g    %9.7g                  !TLATENT HARDENING BY THIS SLIP MODE ON TWIN1, TWIN2:  C IN EQ. (3.28)                         \n',par(oset+5),par(oset+6));
fprintf(fileID,' 0.  0.         				   !TLATENT1 HARDENING BY THIS SLIP MODE ON TWIN1, TWIN2                                          \n');
fprintf(fileID,' %5.1f %5.1f   %5.1f  	               !FOR HPFAC COEF FOR THIS SLIP MODE FOR GRAIN BOUNDARY, TWIN1 BOUNDARY, TWIN2 BOUNDARY          \n',par(oset+1),par(oset+1),par(oset+1));
fprintf(fileID,' 2 5.0  0.0 1.0                    !IQ (LAW FOR Q), Q0,Q1,Q2 FOR LAW FOR Q                                                        \n');
fprintf(fileID,' 2 0.0 0.0 0.0 44.0e+03            !ISHMOD (law for shear modulus),SHEAR MODULUS CONSTANTS                                        \n');oset=15;
fprintf(fileID,'BASAL                                                                                                                             \n');
fprintf(fileID,' 2.95e-10 0.53779 %7.5g            !BURG (m), KOVB3 (J/Km^3), NORM ACTENER g IN EQ. (3.12)                                        \n',par(oset+3));
fprintf(fileID,' %7.5g %7.5g                 !KGENER-K1 IN EQ. (3.8) (1/m), DRAG STRESS-D IN EQ. (3.12) (Pa)                                      \n',par(oset+2),par(oset+4));
fprintf(fileID,' 1.E+07                            !EDOT_O IN EQ. (3.12)                                                                          \n');
fprintf(fileID,' %9.7g 1.E-01                   !INITIAL RHO_S (1/m^2), INITIAL RHO_DEB FOR EACH SLIP MODE (1/m^2)                             \n',par(oset+7));
fprintf(fileID,' 0 %5.1f   0.0   1.0   0.  0.0  0.0 0.0 0.0   !ISL (law for slip),A,B,C,D,E,... SLIP CONSTANTS                                    \n',par(oset));
fprintf(fileID,' 2 %9.7g   %9.7g                     !TLATENT HARDENING BY THIS SLIP MODE ON TWIN1, TWIN2:  C IN EQ. (3.28                          \n',par(oset+5),par(oset+6));
fprintf(fileID,' 0.   0.                           !TLATENT1 HARDENING BY THIS SLIP MODE ON TWIN1, TWIN2                                          \n');
fprintf(fileID,' %5.1f  %5.1f  %5.1f                     !FOR HPFAC COEF FOR THIS SLIP MODE FOR GRAIN BOUNDARY, TWIN1 BOUNDARY, TWIN2 BOUNDARY    \n',par(oset+1),par(oset+1),par(oset+1));
fprintf(fileID,' 2 8.0 0. 1.                     !IQ (LAW FOR Q), Q0,Q1,Q2 FOR LAW FOR Q                                                          \n');
fprintf(fileID,' 2 0.0 0.0 0.0 44.0e+03         !ISHMOD (law for shear modulus),SHEAR MODULUS CONSTANTS                                           \n');
fprintf(fileID,'TWIN1=TENSILE TWINNING                                                                                                            \n');
fprintf(fileID,' 3.0174e-011                       !TWIN BURGERS VECTOR (m)                                                                       \n');
fprintf(fileID,' 1 1                                                                                                                              \n');
fprintf(fileID,' %5.1f 0.0 10.0                                                                                                                   \n',par(33));
fprintf(fileID,' %5.1f 0.0 10.0                                                                                                                   \n',par(33));
fprintf(fileID,' %5.1f                             !HP used in hallpetch term same as just first on in                                            \n',par(34));
fprintf(fileID,'TWIN3=COMPRESSIVE TWINNING                                                                                                        \n');
fprintf(fileID,' 2.725e-011                        !TWIN BURGERS VECTOR (m)                                                                       \n');
fprintf(fileID,' 1 1                                                                                                                              \n');
fprintf(fileID,' %5.1f 0.0 10.0               !B,V,A,mi,sigma                                                                                     \n',par(38));
fprintf(fileID,' %5.1f 0.0 10.0               !B,V,A,mi,sigma                                                                                     \n',par(38));
fprintf(fileID,' %5.1f                             !HP  used in hallpetch term                                                                    \n',par(39)); %Needs to be fit after we get initial yield.
fprintf(fileID,'INTERACTION MATRIX                                                                                                                \n');
fprintf(fileID,' 0.9 0.0 0.0                                                                                                                      \n');
fprintf(fileID,' 0.0 0.9 0.0                                                                                                                      \n');
fprintf(fileID,' 0.0 0.0 0.9                                                                                                                      \n');
fprintf(fileID,'INTERACTION MATRIX                                                                                                                \n');
fprintf(fileID,' 0.0 0.0 0.0                                                                                                                      \n');
fprintf(fileID,' 0.0 0.0 0.0                                                                                                                      \n');
fprintf(fileID,' 0.0 0.0 0.0                                                                                                                      \n');
fprintf(fileID,'INTERACTION MATRIX                                                                                                                \n');
fprintf(fileID,' 0.0 0.0 0.0                                                                                                                      \n');
fprintf(fileID,' 0.0 0.0 0.0                                                                                                                      \n');
fprintf(fileID,' 0.0 0.0 0.0                                                                                                                      \n');
fprintf(fileID,'BARF POWER                                                                                                                        \n');
fprintf(fileID,' %5.1f                                                                                                                             \n',par(4));
fprintf(fileID,'FRACTION OF PARENT WHEN SEC TWINS START CREATING                                                                                  \n');
fprintf(fileID,' %5.1f !TT                                                                                                                          \n',par(37));
fprintf(fileID,' %5.1f !CT                                                                                                                          \n',par(42));
fprintf(fileID,'IPTS TWINTHRES                                                                                                                    \n');
fprintf(fileID,' 0.01                                                                                                                             \n');
fprintf(fileID,'                                                                                                                                  \n');
fprintf(fileID,'-------------------                                                                                                               \n');
fprintf(fileID,'NONSCHMIDT EFFECTS                                                                                                                \n');
fprintf(fileID,'0  7  0  1.0e15  700.0 !flag, nmodes, isch, xepsnsch, tempicnsch                                                                  \n');
fprintf(fileID,' 0.0 0.0 0.0 0.0                                                                                                                  \n');
fprintf(fileID,' 0.0 0.0 0.0 0.0                                                                                                                  \n');
fprintf(fileID,' 0.0 0.0 0.0 0.0                                                                                                                  \n');
fprintf(fileID,' 0.0 0.0 0.0 0.0                                                                                                                  \n');
fprintf(fileID,' 0.0 0.0 0.0 0.0                                                                                                                  \n');
fprintf(fileID,' 0.0 0.0 0.0 0.0                                                                                                                  \n');
fprintf(fileID,' 0.0 0.0 0.0 0.0                                                                                                                  \n');
fprintf(fileID,'                                                                                                                                  \n');
fprintf(fileID,'-------------------                                                                                                               \n');
fprintf(fileID,'SATURATION HARDENING                                                                                                              \n');
fprintf(fileID,'10.0 420.0 220.00 4.25 !tauinit,h0,taus,a                                                                                         \n');
fprintf(fileID,'                                                                                                                                  \n');
fprintf(fileID,'-------------------                                                                                                               \n');
fprintf(fileID,'WORK TO HEAT                                                                                                                      \n');
fprintf(fileID,'16640.0 0.1455e3 0.09544e-1 -68.9e3 0.0 !density[kg/m^3],A0,A1,A2,(C=A0+A1*T+A2/T**2)[J/kg/K],worktoheat                          \n');
fclose(fileID);
end
