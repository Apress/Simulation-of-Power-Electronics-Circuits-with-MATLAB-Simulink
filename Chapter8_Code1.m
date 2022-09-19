clc
clear all
 
syms R1 R2 R D IL VC rC rL VD vIN
 
eq1=-D*R1*IL-(1-D)*R2*IL-R/(R+rC)*VC-(1-D)*VD+D*vIN;
eq2=R/(R+rC)*IL-1/(R+rC)*VC;
 
DC_operatingPoint=solve(eq1,eq2,[IL VC]);
 
disp('IL=')
pretty(simplify(DC_operatingPoint.IL))
 
disp('VC=')
pretty(simplify(DC_operatingPoint.VC))
