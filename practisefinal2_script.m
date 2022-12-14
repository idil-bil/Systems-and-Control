%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% MATLAB PRACTISE FINAL 2 - ELEC 341 - IDIL BIL
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Clean up workspace
clear all; 
clc;

% Student Number
SN    = 21344189;
A = 12;
B = 11;
C = 13;
D = 14;
E = 14;
F = 11;
G = 18;
H = 19;

s = tf('s');            %define s for laplace

%Motor specifications
Rw = A*20*10^(-3);      %(Ohms)
Lw = (B+C)*10^(-3);     %(H)
Jm = F/3;               %(Nms^2/rad)
Bm = H*20;              %(Nms/rad)
Km = D*1/10^3;          %(Nm/A)

%Microcontroller specifications
CF = F*10;              %(Hz)

%Arm specifications
Jb = F * 2;             %(Nms^2/rad)
Bs = G + H;             %(Nms/rad)
Ks = D;                 %(Nm/rad)

%General specifications
tolerance = 0.5;        %(rad/s)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%System Identification

Tr = 33*10^(-3);          %(s) first time it gets to the final value before the peak
Tp = 45*10^(-3);          %(s) first time it gets to the highest value
Ts = 73*10^(-3);          %(s) first time it gets %2 close to settle value after the rise
Kdc = 3.6;                %gain = y value of Tsettle
Pos = (4.26-Kdc)/Kdc*100; %(percent) (Tpeak y - Tsettle y)/Tsettle y * 100

os = Pos/100;                                 %this is more accurate than using the y value of Tpeak
zeta = sqrt(log(os)^2/(pi^2 + log(os)^2));    %(pure) zeta = sqrt[(ln(overshoot)^2/(pi^2+ln(overshoot)^2)]
beta = sqrt(1-zeta^2);                        %(pure) beta = sqrt(1-zeta^2)
Wn = pi/(beta*Tp);                            %(rad/s) omega with Tp = pi/(beta*Tp)
Q1_Ta = minreal(Kdc*(Wn^2)/(s^2+2*zeta*Wn*s+Wn^2));    %(V/V) xfer = Kdc*omega^2/(s^2+2*zeta*omega*s+omega^2)

%state space matrixes for the mechanical system (figure 2)
%F = m*a = B*v = K*x, T = J*sw = B*w = K*q 
%for D -> output number = number of rows, input number = number of columns
Q2_A = [0 0 1 0;
            0 0 0 1;
            -Ks/Jb Ks/Jb -Bs/Jb Bs/Jb;
            Ks/Jm -Ks/Jm Bs/Jm -(Bs+Bm)/Jm];
Q2_B = [0; 0; 0; 1/Jm];
Q2_C = [0 0 0 1;
            1 -1 0 0];
Q2_D = [0; 0];

identity = eye(4);                  %number of elements in the state vector
STM = inv(s*identity-Q2_A);       %system transition matrix = (s*identity matrix - A)^(-1)
G = Q2_C*STM*Q2_B+Q2_D;       %transfer function = C*STM*B+D

%Mechanical admittance (Ym = wm/Tm)
Q3_Ym = G(1);                           %(rad/(Nms))
%Electrical admittance (Ye = Im/Vm)
Q3_Ye = 1/(Rw+s*Lw);                    %(A/V)
%Tf of the motor and arm (Tmech = wm/Vm)
Q3_Tmech = feedback(Q3_Ye*Km*Q3_Ym,Km); %feedback(G,H)

%Ta(V/V) used to convert control voltage into motor voltage
Q4_G = minreal(Q1_Ta*Q3_Tmech*1/s,tolerance);

%tf for the micro-controller daq (H)
Dfb = (2*CF)/(s+2*CF);  %Dfeedback = CF/(Nhat*s+CF) (nhat is not spesifically mentioned so its 1), 2*CF for past years classes
Kfb = 1/1;              %Kfeedback = 1/dc_gain_sensor (there is no sensor so 1)
daq_xfer = Dfb*Kfb;     %data acquisition xfer function is the feedback

Q5_GH = minreal(Q4_G*daq_xfer);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%Design

%controller pole is at 2*CF, as well as the derrivative pole
%use weighted sum to get derivative pole close to controller pole
ws_N_hat = 1;    %weighted sum Nhat = % 2CF/CF = 2 ~1.94 from the table
ws_N = 6;        %weighted sum N -> %corresponding N value to N_hat

%dynamics of controller associated with poles
pid_Dp = (1/s)*(2*CF)/(ws_N_hat*s+2*CF);    %Dpole = 1/s * (2*CF)/(Nhat*s + 2*CF)

%compute Kref (kappa) with Dp
Q6_kappa = margin(pid_Dp*Q5_GH);     %kappa = margin(openloop xfer including Dp)

