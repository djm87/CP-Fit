EPSC Process Input File
Compression in Rolling Direction Displacement Control:
__________________________________
200    Number of steps in process, "nsteps"
1      Process Control Variable, "i_control_var"  [0=temp, 1-6=corresponding known etbc(1-6) or stbc(1-6) component
1      Relative or absolute boundary conditions, i_bc_mode" [0 for relative BC, 1 for absolute BC]
__________________________________
STRAINS
Boundary Conditions of Deformation Step, "ietbc33"
1 1 1
1 0 1
1 1 0
Total Deformation Tensor for the Process, "etbc33"
0.13   0.0   0.0
0.0   0.0   0.0
0.0   0.0   0.0
__________________________________
STRESSES
Boundary Conditions on Stress, "istbc"
0 0 0
  1 0
    1
Total Stress Tensor for the Process, "stbc"
  0.0    0.0   0.0
         0.0   0.0
               0.0
__________________________________
TEMPERATURE
Starting Temperature, "temp_s"			
298.0
Temperature Increment, "deltemp" (change in temperature per step - positive or negative)
0.0
Enforced Temperature Dependence on Elastic Constants (1=Zirconium or 0=Not Zirconium)
0
Define reference values for Macro and Grain Strain Prior to Beginning Process (1=YES, 0=NO)
0              "i_ref_et"
Define reference values for Grain Stress Prior to Beginning Process (1=YES, 0=NO)
0              "i_ref_st"
__________________________________
STRAIN RATE (for dislocation density based haredening law)
0.0001