%This program calculates the small signal transfer 
%functions for Buck converter
R=5;
 
VIN=50;
rin=.1;
 
L=400e-6;
rL=.1;
 
C=100e-6;
rC=.05;
 
rD=.01;
VD=.7;
 
rds=.1;
 
D=.41;
 
R1=rin+rds+rL+R*rC/(R+rC);
R2=rD+rL+R*rC/(R+rC);
 
IL=(R+rC)*(D*VIN-(1-D)*VD)/((R+rC)*R2+R^2+D*(R+rC)*(R1-R2));
 
A=[(R2*(D-1)-R1*D)/L -R/(R+rC)/L;R/(R+rC)/C -1/(R+rC)/C];
B=[(VIN+VD+(R2-R1)*IL)/L D/L;0 0];
CC=[R*rC/(rC+R) R/(R+rC)]; %C shows the capacitance so CC is used for matrix
H=tf(ss(A,B,CC,0));
vO_d=H(1)% transfer function between output voltage and duty ratio
vO_vin=H(2) %transfer function between output voltage and input source
figure(1)
bode(vO_d), grid on
figure(2)
bode(vO_vin), grid on
