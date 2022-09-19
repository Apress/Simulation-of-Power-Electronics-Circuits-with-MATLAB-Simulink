%This program calculates the input and output impedance of the Boost
%converter.
 
clc
 
clear all
syms vg rg d rL L rC C R vC iL rds rD vD io
 
%Converter Dynamical equations
%M1: diL/dt for closed MOSFET.
%M2: dvC/dt for closed MOSFET.
%M3: current of input DC source for closed MOSFET.
%M4: output voltage of converter for closed MOSFET.
 
%M5: diL/dt for open MOSFET.
%M6: dvC/dt for open MOSFET.
%M7: current of input DC source for open MOSFET.
%M8: output voltage of converter for open MOSFET.
 
M1=(-(rg+rL+rds)*iL+vg)/L;
M2=(-vC/(R+rC)+R/(R+rC)*io)/C;
M3=iL;
M4=R/(R+rC)*vC+R*rC/(R+rC)*io;
 
M5=(-(rg+rL+rD+R*rC/(R+rC))*iL-R/(R+rC)*vC-R*rC/(R+rC)*io+vg-vD)/L;
M6=((R/(R+rC))*iL-vC/(R+rC)+R/(R+rC)*io)/C;
M7=iL;
M8=R*rC/(R+rC)*iL-R/(R+rC)*vC+R*rC/(R+rC)*io;
 
%Averaged Equations
diL_dt_ave=simplify(M1*d+M5*(1-d));
dvC_dt_ave=simplify(M2*d+M6*(1-d));
ig_ave=simplify(M3*d+M7*(1-d));
vo_ave=simplify(M4*d+M8*(1-d));
 
%DC Operating Point
DC=solve(diL_dt_ave==0,dvC_dt_ave==0,iL,vC);
IL=DC.iL;
VC=DC.vC;
 
%Linearization
A11=simplify(subs(diff(diL_dt_ave,iL),[iL vC io],[IL VC 0]));
A12=simplify(subs(diff(diL_dt_ave,vC),[iL vC io],[IL VC 0]));
A21=simplify(subs(diff(dvC_dt_ave,iL),[iL vC io],[IL VC 0]));
A22=simplify(subs(diff(dvC_dt_ave,vC),[iL vC io],[IL VC 0]));
AA=[A11 A12;A21 A22];
 
B11=simplify(subs(diff(diL_dt_ave,io),[iL vC io],[IL VC 0]));
B12=simplify(subs(diff(diL_dt_ave,vg),[iL vC io],[IL VC 0]));
B13=simplify(subs(diff(diL_dt_ave,d),[iL vC io],[IL VC 0]));
 
B21=simplify(subs(diff(dvC_dt_ave,io),[iL vC io],[IL VC 0]));
B22=simplify(subs(diff(dvC_dt_ave,vg),[iL vC io],[IL VC 0]));
B23=simplify(subs(diff(dvC_dt_ave,d),[iL vC io],[IL VC 0]));
 
BB=[B11 B12 B13;B21 B22 B23];
 
C11=simplify(subs(diff(ig_ave,iL),[iL vC io],[IL VC 0]));
C12=simplify(subs(diff(ig_ave,vC),[iL vC io],[IL VC 0]));
 
C21=simplify(subs(diff(vo_ave,iL),[iL vC io],[IL VC 0]));
C22=simplify(subs(diff(vo_ave,vC),[iL vC io],[IL VC 0]));
CC=[C11 C12; C21 C22];
 
D11=simplify(subs(diff(ig_ave,io),[iL vC io],[IL VC 0 ]));
D12=simplify(subs(diff(ig_ave,vg),[iL vC io],[IL VC 0]));
D13=simplify(subs(diff(ig_ave,d),[iL vC io],[IL VC 0]));
 
D21=simplify(subs(diff(vo_ave,io),[iL vC io],[IL VC 0 ]));
D22=simplify(subs(diff(vo_ave,vg),[iL vC io],[IL VC 0]));
D23=simplify(subs(diff(vo_ave,d),[iL vC io],[IL VC 0]));
DD=[D11 D12 D13;D21 D22 D23];
 
%Components Values
%Variables have underline are used to store the numeric values of components
%Variables without underline are symbolic variables.
%for example:
%L: symbolic vvariable shows the inductor inductance
%L_: numeric variable  shows the inductor inductance value.
L_=120e-6;
rL_=.01;
C_=100e-6;
rC_=.05;
rds_=.04;
rD_=.01;
VD_=.7;
D_=.6;
VG_=12;
rg_=.1;
R_=50;
 
AA_=eval(subs(AA,[vg rg rds rD vD rL L rC C R d io],[VG_ rg_ rds_ rD_ VD_ rL_ L_ rC_ C_ R_ D_ 0]));
BB_=eval(subs(BB,[vg rg rds rD vD rL L rC C R d io],[VG_ rg_ rds_ rD_ VD_ rL_ L_ rC_ C_ R_ D_ 0]));
CC_=eval(subs(CC,[vg rg rds rD vD rL L rC C R d io],[VG_ rg_ rds_ rD_ VD_ rL_ L_ rC_ C_ R_ D_ 0]));
DD_=eval(subs(DD,[vg rg rds rD vD rL L rC C R d io],[VG_ rg_ rds_ rD_ VD_ rL_ L_ rC_ C_ R_ D_ 0]));
 
sys=ss(AA_,BB_,CC_,DD_);
sys.stateName={'iL','vC'};
sys.inputname={'io','vg','d'};
sys.outputname={'ig','vo'};
 
ig_io=sys(1,1);
ig_vg=sys(1,2);
ig_d=sys(1,3);
 
vo_io=sys(2,1);
vo_vg=sys(2,2);
vo_d=sys(2,3);
 
Zin=1/ig_vg; %input impedance
Zout=vo_io;  %output impedance 
 
%Draws the bode diagram of input/output impedance
figure(1)
bode(Zin), grid minor
 
figure(2)
bode(Zout), grid minor
 
%Display the DC operating point of converter
disp('steady state operating point of converter')
disp('IL')
disp(eval(subs(IL,[vg rg rds rD vD rL L rC C R d io],[VG_ rg_ rds_ rD_ VD_ rL_ L_ rC_ C_ R_ D_ 0])));
disp('VC')
disp(eval(subs(VC,[vg rg rds rD vD rL L rC C R d io],[VG_ rg_ rds_ rD_ VD_ rL_ L_ rC_ C_ R_ D_ 0])));
