%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% MATLAB PRACTISE FINAL 1 - ELEC 341 - IDIL BIL
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

s = tf('s');                    %define s for laplace

%Motor specifications
Vin = 12;                       %(V)
Io = (250 + (10*A))*10^(-3);    %(A) no load current (Io = Iw if there is no load)
Rw = 0.2;                       %(Ohms)
Lw = 0.5*10^(-3);               %(H)
Jr = (100 + B)*10^(-7);         %(kgm^2)
Km = 1/[(600 + (10*C))*pi/30];  %(RPM/V to rad/(sV)), unit make more sence when it has 1/A or 1/V

%Microcontroller specifications
CF = 100 + (10*A);              %(Hz)

%Load specifications
Jl = (A + B + C)*10^(-6);       %(Nms^2/rad)
Bl = (D + E + F)*10^(-6);       %(Nms^2/rad)
Kl1 = 2*10^(-3);                %(Nms^2/rad)
Kl2 = 3 * G * 10^(-3);          %(Nms^2/rad)

%Amplifier specifications
amp_DC_gain = B;                %(V/V)
amp_wn = A * B * C * D;         %(rad/s)
amp_zeta = D/10;                %(pure)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%System Identification

%tf for the voltage amplifier (part of G)
amp_xfer = amp_DC_gain*(amp_wn^2)/(s^2+2*amp_zeta*amp_wn*s+amp_wn^2);   %(V/V) xfer = Kdc*omega^2/(s^2+2*zeta*omega*s+omega^2)

%tf for the micro-controller daq (H)
Dfb = (2*CF)/(s+2*CF);  %Dfeedback = CF/(Nhat*s+CF) (nhat is not spesifically mentioned so its 1), 2*CF for past years classes
Kfb = 1/1;              %Kfeedback = 1/dc_gain_sensor (there is no sensor so 1)
daq_xfer = Dfb*Kfb;     %data acquisition xfer function is the feedback

%to find Br
wo = (Vin-Rw*Io)/Km;
Br = Km*Io/wo;

%state space matrixes
matr_A = [-Rw/Lw -Km/Lw 0;
            Km/(Jr+Jl) -(Bl+Br)/(Jr+Jl) -(Kl1*Kl2/(Kl1+Kl2))/(Jr+Jl);
            0 1 0];
matr_B = [1/Lw; 0; 0];
matr_C = [0 0 1];
matr_D = [0];
%for D -> output number = number of rows, input number = number of columns

identity = eye(3);                  %number of elements in the state vector
STM = inv(s*identity-matr_A);       %system transition matrix = (s*identity matrix - A)^(-1)
q_xfer = matr_C*STM*matr_B+matr_D;  %(part of G) transfer function for angle = C*STM*B+D

GH = amp_xfer*daq_xfer*q_xfer;       %GH

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%Design

%controller pole is at 2*CF, as well as the derrivative pole
%use weighted sum to get derivative pole close to controller pole
ws_N_hat = 1;    %weighted sum Nhat = % 2CF/CF = 2 ~1.94 from the table
ws_N = 6;        %weighted sum N -> %corresponding N value to N_hat

%dynamics of controller associated with poles
pid_Dp = (1/s)*(2*CF)/(ws_N_hat*s+2*CF);    %Dpole = 1/s * (2*CF)/(Nhat*s + 2*CF)

%compute Kref (kappa) with Dp
Kref = margin(pid_Dp*GH);     %kappa = margin(openloop xfer including Dp)

%find wxo (frequency where we start trying numbers to find the max phase margin) usign Kref
[G_m P_m wxo_g wxo_p] = margin(Kref*pid_Dp*GH);   %gain margin, phase margin, freq gain, freq phase (freq gain = freq phase when margin is done with Kref)

%Search for initial zero positions
if 0    %starts from the inital frequency found and finds the real zeros that give the max phase margin
    [pid_D_rec_r2, zero12_r, maxPm_r2] = maxPmRealZeros(wxo_g,Kref,pid_Dp,GH);    %(initial freq, Kref kappa, Dp, GH after Dp)
end     %[D, zero matrix[zero1, zero2], max phase margin]
    
if 0    %starts from the inital frequency found and finds the complex zeros that give the max phase margin
    [pid_D_rec_i1, zero12_i1, maxPm_i1, img_search_valid1, invalid_matrix1] = maxPmComplexZeros(wxo_g,Kref,pid_Dp,GH);    %(initial freq, Kref kappa, Dp, GH after Dp)
end     %[D, zero matrix[zero1, zero2], max phase margin, if 1: the code found the max phase margin if 0: decrease resolution to get 1, gives the invalid matrix so change resolution or iteration number]
%if it doesn't matter if the zeros are real or complex compare the phase margins in the terminal and choose the bigger ones zeros
    
%hard code these values after running so you don't need to run again
% zero12_r = [-5.6346 -5.6356];
zero12_i1 = [-3.9721 + 3.8765i -3.9721 - 3.8765i];
% maxPm_r2 = 99.4645;
maxPm_i1 = 109.6686;
%choose the one that gives the bigger phase margin (Pm)
% pid_D_rec_r2 = (13.86*s^2 + 156.2*s + 440)/(s^2 + 440*s);
pid_D_rec_i1 = (14.28*s^2 + 113.5*s + 440)/(s^2 + 440*s);
pid_D_maxPm = minreal(pid_D_rec_i1);  %complex in this case

%poles and zeros for initial designed controller
pid_pole = -2*CF;           %always -2*CF
pid_zero1 = zero12_i1(1);   %from the chosen zeros
pid_zero2 = zero12_i1(2);   %from the chosen zeros

%dynamics gains
init_pid_Ki = 1;                                                         %Ki = 1 (for initial)
init_pid_Kp = 1/pid_pole-(pid_zero1+pid_zero2)/(pid_zero1*pid_zero2);    %Kp = 1/p - (z1+z2)/(z1*z2)
init_pid_Kd = 1/(pid_zero1*pid_zero2)+init_pid_Kp/pid_pole;              %Kd when Kp is defined = 1/(z1*z2) + Kp/p

%desired hase margin, if its not given define as 45
desiredPm = 30;

%master gain of the controller for a desired phase argin
init_pid_masterK = getK4PhaseMargin(Kref,pid_D_maxPm,GH,desiredPm);     %master gain (K) = (Kref, D, openloop after D, desired phase margin)

%%%
% REQUIREMENTS
% OS > 2 %
%%%
% CONSTRAINTS
% do not change CF
% do not change FDD filter
%%%
% GOALS
% No undershoot
% OS small
% Ts small
%%%

%tune the following gains to satisfy RCGs listed above
Ki_factor = 0.65;  
tuned_pid_Ki = Ki_factor;
tuned_pid_Kp = init_pid_Kp * 0.95;
tuned_pid_Kd = init_pid_Kd * 0.7;
tuned_pid_masterK = init_pid_masterK * 0.70;

%full dynamics of the controller 
tuned_pid_D = tuned_pid_Kp + tuned_pid_Ki*(1/s) + tuned_pid_Kd*(-pid_pole*s)/(1*s-pid_pole);

%the response with the tuned controller
tuned_response = feedback(tuned_pid_masterK*tuned_pid_D*amp_xfer*q_xfer,daq_xfer);  %feedback(G, H)
[gain_m, phase_m] = margin(tuned_pid_masterK*tuned_pid_D*amp_xfer*q_xfer,daq_xfer); %gain and phase margin of the tuned version

%plot the response
figure(100); clf;
step(tuned_response);
grid on;
stepinfo(tuned_response)