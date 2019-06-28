function [ output_args ] = WriteSxFile(fname,systemPar,dir)
    % Parameters to write 
increments=systemPar{8}(2)*1000;
fileID=fopen(fname,'w');
if dir=='C1'
    fprintf(fileID,'%d   1   0.001   298.         nsteps  ictrl  eqincr  temp\n',increments);
    fprintf(fileID,'* boundary conditions\n');
    fprintf(fileID,'    1       1       1           iudot    |    flag for vel.grad.\n');
    fprintf(fileID,'    1       0       1                    |    (0:unknown-1:known)\n');
    fprintf(fileID,'    1       1       0                    |\n');
    fprintf(fileID,'                                         |\n');
    fprintf(fileID,'   -0.001    0.      0.           udot   |    vel.grad\n');
    fprintf(fileID,'    0.    0.0005   0.                    |\n');
    fprintf(fileID,'    0.      0.    0.0005                 |\n');
    fprintf(fileID,'                                         |\n');
    fprintf(fileID,'    0       0       0           iscau    |    flag for Cauchy\n');
    fprintf(fileID,'            1       0                    |\n');
    fprintf(fileID,'                    1                    |\n');
    fprintf(fileID,'                                         |\n');
    fprintf(fileID,'    0.      0.      0.          scauchy  |    Cauchy stress\n');
    fprintf(fileID,'            0.      0.                   |\n');
    fprintf(fileID,'                    0.                   |\n');
elseif dir == 'C2'
    fprintf(fileID,'%d   2   0.001   298.         nsteps  ictrl  eqincr  temp\n',increments);
    fprintf(fileID,'* boundary conditions\n');
    fprintf(fileID,'    0       1       1           iudot    |    flag for vel.grad.\n');
    fprintf(fileID,'    1       1       1                    |    (0:unknown-1:known)\n');
    fprintf(fileID,'    1       1       0                    |\n');
    fprintf(fileID,'                                         |\n');
    fprintf(fileID,'   0.0005    0.      0.           udot   |    vel.grad\n');
    fprintf(fileID,'    0.    -0.001   0.                    |\n');
    fprintf(fileID,'    0.      0.    0.0005                 |\n');
    fprintf(fileID,'                                         |\n');
    fprintf(fileID,'    1       0       0           iscau    |    flag for Cauchy\n');
    fprintf(fileID,'            0       0                    |\n');
    fprintf(fileID,'                    1                    |\n');
    fprintf(fileID,'                                         |\n');
    fprintf(fileID,'    0.      0.      0.          scauchy  |    Cauchy stress\n');
    fprintf(fileID,'            0.      0.                   |\n');
    fprintf(fileID,'                    0.                   |\n');
elseif dir == 'C3' 
    fprintf(fileID,'%d   3   0.001   298.         nsteps  ictrl  eqincr  temp\n',increments);
    fprintf(fileID,'* boundary conditions\n');
    fprintf(fileID,'    0       1       1           iudot    |    flag for vel.grad.\n');
    fprintf(fileID,'    1       0       1                    |    (0:unknown-1:known)\n');
    fprintf(fileID,'    1       1       1                    |\n');
    fprintf(fileID,'                                         |\n');
    fprintf(fileID,'   0.0005    0.      0.           udot   |    vel.grad\n');
    fprintf(fileID,'    0.    0.0005   0.                    |\n');
    fprintf(fileID,'    0.      0.    -0.001                 |\n');
    fprintf(fileID,'                                         |\n');
    fprintf(fileID,'    1       0       0           iscau    |    flag for Cauchy\n');
    fprintf(fileID,'            1       0                    |\n');
    fprintf(fileID,'                    0                    |\n');
    fprintf(fileID,'                                         |\n');
    fprintf(fileID,'    0.      0.      0.          scauchy  |    Cauchy stress\n');
    fprintf(fileID,'            0.      0.                   |\n');
    fprintf(fileID,'                    0.                   |\n');
end
fclose(fileID);
end
