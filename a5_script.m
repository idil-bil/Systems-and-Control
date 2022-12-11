%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% MATLAB ASSIGNMENT 5 - ELEC 341 - IDIL BIL
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

%Defining Variables of the Mechanical Circuit from assignment 3
%F = m*a = B*v = K*x, T = J*sw = B*w = K*theta 
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%Q1: State Matrixes

%[sd; sv] = A * [d; v] + B * [F0; ](inputs)
Q1.A = [0 0 0 0 1 0 0 0; 
        0 0 0 0 0 1 0 0; 
        0 0 0 0 0 0 1 0; 
        0 0 0 0 0 0 0 1;
        -K20/M0 0 K20/M0 0 (-B20-B0)/M0 0 B20/M0 0;
        0 (-K21-K1)/M1 K21/M1 0 0 -B31/M1 0 B31/M1;
        K20/M2 K21/M2 (-K32-K20-K21)/M2 K32/M2 B20/M2 0 -B20/M2 0;
        0 0 K32/M3 -K32/M3 0 B31/M3 0 -B31/M3];
Q1.B = [0; 0; 0; 0; 300/M0; 0; 0; 0];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%Q2: Output Matrices

%[d3; fK20](outputs) = C * [d; v] + D * [F0; ](inputs)    
Q2.C = [0 0 0 1 0 0 0 0; K20 0 -K20 0 0 0 0 0];
Q2.D = [0; 0];

%plots to check Q1 and Q2
if 0
        sys = ss(Q1.A,Q1.B,Q2.C,Q2.D);
        impulse(sys);   %compare with plots from assignment 3
        grid on;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%Q3: Actuated System

%Variable Declaration
Rw = A/3;        %(ohms)
Lw = B*10^(-3);  %(henrys)
Km = G/2;        %(N/a or Vs/m)
Vin = 150;       %(volts)
M0 = A+B;        %(kg)
B0 = C+D;        %(Ns/m)
N0 = 1:2;        %speed ratio

%[sd; sv; sIw] = A * [d; v; Iw] + B * [Vin; ](input)
A = [0 0 0 0 1 0 0 0 0; 
     0 0 0 0 0 1 0 0 0; 
     0 0 0 0 0 0 1 0 0; 
     0 0 0 0 0 0 0 1 0;
     -K20/M0 0 K20/M0 0 (-B20-B0)/M0 0 B20/M0 0 Km/(2*M0);
     0 (-K21-K1)/M1 K21/M1 0 0 -B31/M1 0 B31/M1 0;
     K20/M2 K21/M2 (-K32-K20-K21)/M2 K32/M2 B20/M2 0 -B20/M2 0 0;
     0 0 K32/M3 -K32/M3 0 B31/M3 0 -B31/M3 0;
     0 0 0 0 -Km/Lw 0 0 0 -Rw/Lw];
B = [0; 0; 0; 0; 0; 0; 0; 0; 1/Lw];

%[d3; ](outputs) = C * [d; v; Iw] + D * [Vin; ](input)
C = [0 0 0 1 0 0 0 0 0];
D = [0];
 
identity = eye(9);
s = tf('s');              %define s for laplace
STM = inv(s*identity-A);  %system transition matrix = (s*identity matrix - A)^(-1)
Q3.G = C*STM*B+D;         %transfer function = C*STM*B+D

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%Q4: Gravity
%gravity only affects inputs to B and D changes

%[sd; sv; sIw] = A * [d; v; Iw] + B * [Vin; Gravity](input)
B4 = [0 0; 0 0; 0 0; 0 0; 0 -9.81; 0 -9.81; 0 -9.81; 0 -9.81; 1/Lw 0];

%[d3; ](outputs) = C * [d; v; Iw] + D * [Vin; Gravity](input)
%for D -> output number = number of rows, input number = number of columns
D4 = [0 0];

STM4 = inv(s*identity-A);   %system transition matrix = (s*identity matrix - A)^(-1)
TF4 = C*STM4*B4+D4;         %transfer function = C*STM*B+D
Q4.Gv = TF4(1);
Q4.Gg = TF4(2);
