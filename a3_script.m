%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% MATLAB ASSIGNMENT 3 - ELEC 341 - IDIL BIL
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Clean up workspace
clear all; 
clc;

% Student Number
SN    = ABCDEFGH;
A = 12;
B = 11;
C = 13;
D = 14;
E = 14;
F = 11;
G = 18;
H = 19;

%Defining Variables of the Mechanical Circuits
%F = ma = Bv = Kd, T = Jw = Bw = K*theta 
M0 = A;   %(kg)
M3 = A;   %(kg)
M1 = B;   %(kg)
M2 = B;   %(kg)
B0 = C;   %(Ns/m)
B1 = D;   %(Ns/m)
B20 = E;  %(Ns/m)
B31 = F;  %(Ns/m)
K1 = G;   %(N/m)
K32 = G;  %(N/m)
K20 = H;  %(N/m)
K21 = H;  %(N/m)
F0 = 300; %(N)

%Define Variables for the Electrical Circuits
%R = 1/B, L = 1/K, C = M, I = F, V = v
C0 = M0;     %(farads)
C3 = M3;     %(farads)
C1 = M1;     %(farads)
C2 = M2;     %(farads)
R0 = 1/B0;   %(ohms)
R1 = 1/B1;   %(ohms)
R20 = 1/B20; %(ohms)
R31 = 1/B31; %(ohms)
L1 = 1/K1;   %(henrys)
L32 = 1/K32; %(henrys)
L20 = 1/K20; %(henrys)
L21 = 1/K21; %(henrys)
I0 = F0;     %(amperes)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%Q1: Mass Distance

%Admittances
s = tf('s');    %define s for laplace
yC0 = s*C0;
yC3 = s*C3;
yC1 = s*C1;
yC2 = s*C2;
yR0 = 1/R0;
yR1 = 1/R1;
yR20 = 1/R20;
yR31 = 1/R31;
yL1 = 1/(s*L1);
yL32 = 1/(s*L32);
yL20 = 1/(s*L20);
yL21 = 1/(s*L21);

%Admittance Matrix
%written after converting the system and using circuit analysis with matrix method
Y = [yR0+yC0+yL20+yR20 0 -yL20-yR20 0;
    0 yL1+yL21+yR31+yC1 -yL21 -yR31; 
    -yL20-yR20 -yL21 yL32+yL20+yR20+yL21+yC2 -yL32; 
    0 -yR31 -yL32 yL32+yR31+yC3];

%Input Current Matrix
%take input as 1 to find transfer functions
I = [1; 0; 0; 0];   %the current source only affects V0

%Voltage Transfer Function Matrix
V = inv(Y)*I;   %I = Y*V

%Transfer Function for V3
V3_tf = V(4);   %V = [V0; V1; V2; V3]

%Transfer Function for d3
d3_tf = V3_tf*1/s;      %use tf input (1) and integrate velocity in s domain (1/s) to get distance from velocity
Q1.G = minreal(d3_tf);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%Q2: Spring Force

iL20_tf = (V(1)-V(3)) / (s*L20);    %xfer is found with force which is the current for electronic systems
Q2.G = minreal(iL20_tf,0.0001);     %minreal percision

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%Q3: Mass Distance

%Admittance Matrix but spring K1 replaced with damper B1
q3_Y = [yR0+yC0+yL20+yR20 0 -yL20-yR20 0;
    0 yR1+yL21+yR31+yC1 -yL21 -yR31; 
    -yL20-yR20 -yL21 yL32+yL20+yR20+yL21+yC2 -yL32; 
    0 -yR31 -yL32 yL32+yR31+yC3];

%Voltage Transfer Function Matrix but spring K1 replaced with damper B1
q3_V = inv(q3_Y)*I;    %I = inv(Y)*V

%Transfer Function for V3 but spring K1 replaced with damper B1
q3_V3 = q3_V(4);       %V = [V0; V1; V2; V3]

%Transfer Function for d3 but spring K1 replaced with damper B1
q3_d3 = q3_V3*1/s;     %scale by input value (1) and integrate velocity in s domain (1/s) to get distance from velocity
Q3.G = minreal(q3_d3);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%Q4: Rotor Speed

%Motor Model Variables
Rw = 1 + A/10;                %(ohms)
Lw = (100 + 10*B) * 10^(-6);  %(henrys)
Jr = C/10 * 10^(-6);          %(Nms^2)
Br = (D+E+F) * 10^(-6);       %(Nms)
Km = (10+G) * 10^(-3);        %(Vs)
Km_p = (10+G) * 10^(-3);      %(Nm/A) prime
Vin = 12;                     %(volts)

%nodal analysis after both circutis are electrical systems leads to 2 eqns with 2 unknowns Iw and w
%through algebraic methods it is possible to find an expression for w in terms of known variables
q4_w = (Km_p*Vin/(s*Lw+Rw))/((Km_p*Km/(s*Lw+Rw))+Br+s*Jr);
q4_xfer = q4_w/Vin;             %tf = Vout/Vin
Q4.G = minreal(q4_xfer);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%Q5: Winding Current

% using the same equations and substituting w value from q4 gives the following
q5_Iw = (Vin-Km*q4_w)/(s*Lw+Rw);
q5_Iw = q5_Iw/Vin;              %tf = Vout/Vin
Q5.G = minreal(q5_Iw,0.0001);
