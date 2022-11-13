%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% MATLAB ASSIGNMENT 4 - ELEC 341 - IDIL BIL
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

%Defining Variables for the MIMO Block Diagram 
s = tf('s');    %define s for laplace
G1 = 1/(s+A);
G2 = 10/(s+B);
G3 = 10/(s+C);
G4 = 10/(s+D);
G5 = 1/(s+E);
H1 = 100/(s+F);
H2 = 1000/(s+G);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%Q1: Block Diagram Manipulation
%to check these on simulink draw the diagram with the library and check the scope graph
%summation is under math, tf is under continuous, step (input) is under source (step time = 0), scope (output) is under commonly used

Q1.G11 = minreal((G1*G2)/(1+G1*(G5*H2*(G3+G2*H1*G4))));

%alternate way
%syms s Y1 Y2 U1 U2 
%U2 = 0; 
%eq1 = (U1 - Y2*H2)*G1*G2 == Y1; 
%eq2 = ((U1 - Y2*H2)*G1*G3 + (Y1*H1-U2)*G4)*G5 == Y2; 
%[Y1, Y2] = solve(eq1, eq2, Y1, Y2); 
%G11 = Y1/U1; 
%[num11, den11] = numden(G11); 
%tfnum11 = sym2poly(num11); 
%tfden11 = sym2poly(den11); 
%Q1.G11 = minreal(tf(tfnum11,tfden11));

G12 = (G1*G3+G1*G2*H1*G1)*G5;
H12 = H2;
Q1.G12 = minreal(G12/(1+G12*H12));  %G/(1+G*H)

Q1.G21 = minreal(((G4*G2*G5*H2*G1)/(1+G3*G5*H2*G1))/(1+((G4*G2*G5*H2*G1)/(1+G3*G5*H2*G1))*H1));

Q1.G22 = minreal(-G5*G4/(1+G5*(-G1*H2)*(H1*G2*G4+G3)));

%plots to check Q1
if 0
    figure(1); clf;
    step(Q1.G11)    %want the the scope graph to be the same as this
    figure(2); clf;
    step(Q1.G12)    %want the the scope graph to be the same as this
    figure(3); clf;
    step(Q1.G21)    %want the the scope graph to be the same as this
    figure(4); clf;
    step(Q1.G22)    %want the the scope graph to be the same as this
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%Q2: MISO System

%superposition
Q2.G1 = minreal(Q1.G11 + Q1.G21);
Q2.G2 = minreal(Q1.G12 + Q1.G22);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%Q3: Fan Motor

%variable declaration
Jf = G/30*10^(-6);            %(Nms^2)
Br = (D+E+F) * 10^(-6);       %(Nms)
Bf = H*Br*10^(-6);            %(Nms)
Rw = 1 + A/10;                %(ohms)
Lw = (100 + 10*B) * 10^(-6);  %(henrys)
Jr = C/10 * 10^(-6);          %(Nms^2)
Km = (10+G) * 10^(-3);        %(Vs)
Km_p = (10+G) * 10^(-3);       %(Nm/A) prime

s = tf('s');            %define s for laplace 
Jtot = Jr*Jf/(Jr+Jf);   %series like capacitors
Btot = Br*Bf/(Br+Bf);   %series like resistors

Q3.Ye = minreal(1/(s*Lw + Rw));    %admittance of resistance and inductance
Q3.Ym = minreal(1/Btot + s*Jtot);  %admittance of B(resistance) and J(capacitance)

fb = feedback((Ye*Km*Ym), Km);     %feedback(G, H) = G/(1+G*H)
Q3.G = minreal(fb);
