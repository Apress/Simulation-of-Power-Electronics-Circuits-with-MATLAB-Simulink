% This program extracts the small signal transfer function
clc
clear all;
% Elements values
R=5;     %Load resistor                  
VIN=50;  %Input source voltage 
rin=.1;  %Input source internal resistance
L=400e-6;%inductor
rL=.1;   %inductor series resistance
C=100e-6;%capacitor
rC=.05;  %capacitor series resistance
rD=.01;  %Diode series resistance
VD=.7;   %Diode forward voltage drop
rds=.1;  %MOSFET on resistance 
D=.41;   %Duty ratio
 
% Symbolic variables
%iL: inductor current 
%vC: capacitor voltage
%vin: input voltage source
%vD: diode forward voltage drop
%d: duty cycle
syms iL vC vin vD d
 
%CLOSED MOSFET EQUATIONS
M1=(-(rin+rds+rL+(R*rC/(R+rC)))*iL-R/(R+rC)*vC+vin)/L;%d(iL)/dt for closed MOSFET
M2=(R/(R+rC)*iL-1/(R+rC)*vC)/C;                       %d(vC)/dt for closed MOSFET 
vO1=R*(rC/(rC+R)*iL+1/(R+rC)*vC);
%OPENED MOSFET EQUATIONS
M3=(-(rD+rL+R*rC/(R+rC))*iL-R/(R+rC)*vC-vD)/L;        %d(iL)/dt for opened MOSFET
M4=(R/(R+rC)*iL-1/(R+rC)*vC)/C;                       %%d(vC)/dt for opened MOSFET 
vO2=R*(rC/(rC+R)*iL+1/(R+rC)*vC);
%AVERAGING
MA1= simplify(d*M1+(1-d)*M3);
MA2= simplify(d*M2+(1-d)*M4);
vO= simplify(d*vO1+(1-d)*vO2);
% DC OPERATING POINT CALCULATION
MA_DC_1=subs(MA1,[vin vD d],[VIN VD D]);
MA_DC_2=subs(MA2,[vin vD d],[VIN VD D]);
 
DC_SOL= solve(MA_DC_1==0,MA_DC_2==0,iL,vC);
 
IL=eval(DC_SOL.iL);   %IL is the inductor current steady state value 
VC=eval(DC_SOL.vC);   %VC is the capacitor current steady state value
 
%LINEARIZATION
% .
% x=Ax+Bu
%vector x=[iL;vC] is assumed. vector x is states.
%u=[vin;d] where vin=input voltage source and d=duty. vector u is system inputs.
%
A11=subs(simplify(diff(MA1,iL)),[iL vC d vD],[IL VC D VD]);
A12=subs(simplify(diff(MA1,vC)),[iL vC d vD],[IL VC D VD]);
 
A21=subs(simplify(diff(MA2,iL)),[iL vC d vD],[IL VC D VD]);
A22=subs(simplify(diff(MA2,vC)),[iL vC d vD],[IL VC D VD]);
 
A=eval([A11 A12;
        A21 A22]);    %variable A is matrix A in state space equation
    
 
B11=subs(simplify(diff(MA1,vin)),[iL vC d vD vin],[IL VC D VD VIN]);
B12=subs(simplify(diff(MA1,d)),[iL vC d vD vin],[IL VC D VD VIN]);
 
B21=subs(simplify(diff(MA2,vin)),[iL vC d vD vin],[IL VC D VD VIN]);
B22=subs(simplify(diff(MA2,d)),[iL vC d vD vin],[IL VC D VD VIN]);
 
 
B=eval([B11 B12;    
        B21 B22]);    % variable B is matrix B in state space equation
 
CC1=subs(simplify(diff(vO,iL)),[iL vC d vD],[IL VC D VD]);
CC2=subs(simplify(diff(vO,vC)),[iL vC d vD],[IL VC D VD]);
CC=eval([CC1 CC2]);   %variable CC is matrix C in state space equation
                      % variable D shows duty so DD is used. 
DD11=subs(simplify(diff(vO,vin)),[iL vC d vD vin],[IL VC D VD VIN]);
DD12=subs(simplify(diff(vO,d)),[iL vC d vD vin],[IL VC D VD VIN]);
 
DD=eval([DD11 DD12]); % variable DD is matrix D in state space equation
                      % variable D shows duty so DD is used.  
 
H=tf(ss(A,B,CC,DD));
               
               %transfer function between input source and load resistor voltage
               % ~               
vR_vin=H(1,1)  % vR(s)
               % ----                
               % ~     
               % vin(s)     
               
               %transfer function between duty ratio and load resistor voltage
               %~               
vR_d=H(1,2)    %vR(s)
               %----                
               %~     
               %d(s)
