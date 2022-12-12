%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% MATLAB ASSIGNMENT 8 - ELEC 341 - IDIL BIL
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

%System Specs from assignment 6
%pole-zero declarations for EMS (G - controller) and SEN (H - feedback) (figure 2)
ems_p1 = -A;            %X
ems_p2 = -3*B;          %X, most dominant pole for ems
ems_p3 = -2*D+2*D*i;    %X
ems_p4 = -2*D-2*D*i;    %X
ems_z1 = -5*C;          %O
ems_dcgain = 300/G;     %gain when s = 0 (X)
sen_p1 = -25*E;         %X
sen_dcgain = 1;         %gain when s = 0 (x)

%xfer declarations
%zpk([zeros],[poles],gain), (s-zeros)/(s-poles)
ems_xfer = zpk([ems_z1],[ems_p1,ems_p2,ems_p3,ems_p4],1);   %1 is the coefficient not gain (K)
sen_xfer = zpk([],[sen_p1],1);                              %1 is the coefficient not gain (K)

%gain adjustments
%overall dc gain/dynamic dc gain = coefficient in the xfer function 
ems_gain = ems_dcgain/dcgain(ems_xfer);                          %gain = dc gain/coefficient   actual K = X/Y
ems_xfer = zpk([ems_z1],[ems_p1,ems_p2,ems_p3,ems_p4],ems_gain);
sen_gain = sen_dcgain/dcgain(sen_xfer);                          %gain = dc gain/coefficient   actual K = X/Y
sen_xfer = zpk([],[sen_p1],sen_gain);

%feedback parameters (figure 1)
Kf = 1/sen_dcgain;  %Kf = 1/gain of the sensor before (X)
CF = 10*F;          %(Hz)
s = tf("s");        %define s for laplace
Df = CF/(s+CF);     %Dfeedback = CF/[Nhat*s +CF] Nhat = 1 if there are on weighted sum filters (also not used when CF is smaller than 10 times the leftmost pole)

%open loop xfer
ol_xfer = ems_xfer*sen_xfer*Df*Kf;     %K is not inluded since its not calculated yet  (GH)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%Q1: Lag Dynamics

%cancel out the most dominant pole
mdp = ems_p1;
pi_zero = mdp;

%controller parameters
pi_Kp = -1/pi_zero;     %Kp = -1/zero in PI controllers
pi_Ki = 1;              %Ki = 1 in PI controllers

%dynamics of the PI Controller
pi_D = (pi_Ki*1/s+pi_Kp);       %D = Ki*1/s + Kp

Q1.Kp = pi_Kp;
Q1.D = pi_D;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%Q2: Tuned PI Control

%find the ultimate gain Ku(kappa)
pi_Ku = margin(Q1.D*ol_xfer);

if 0
    max_OS_value = 1.1;     %overshoot requirement is 10%
    minimum_Ts = 10000000;

    for cnt1 = -1:-1:-100
        for cnt2 = 1:1:95
            valid = 1;
            testZ = cnt1; %will get incremented
            testK = cnt2; 

            %adjust new xfer for test
            pi_Kp_test = -1/testZ;
            pi_Ki_test = 1;
            pi_D_test = (pi_Ki_test*1/s+pi_Kp_test);

            %find new kappa value
            pi_Ku_test = margin(pi_D_test*ol_xfer);

            %try K values
            pi_K_test = pi_Ku_test*testK/100;

            %get CL xfer
            pi_clxfer_test = feedback(pi_K_test*pi_D_test*ems_xfer,sen_xfer*Kf*Df);

            %check peak value
            [ytest xtest] = step(pi_clxfer_test);

            for cnt3 = 1:length(ytest)
                if(ytest(cnt3) > max_OS_value)
                    valid = 0;
                    break;
                end
            end

            %check min Ts
            if(valid == 1)
                if(stepinfo(pi_clxfer_test).SettlingTime < minimum_Ts)      %change the "SettlingTime" to "Overshoot" if you want to find min overshoot
                    minimum_Ts = stepinfo(pi_clxfer_test).SettlingTime;
                    pi_design_K = pi_K_test;
                    pi_design_Z = testZ;
                end
            end
        end
    end
end

%hard coded design values
pi_design_K = 0.4706;                   %get this value from the loop
pi_design_Z = -10;                      %get this value from the loop
pi_design_D = (1*1/s-1/pi_design_Z);    %get this value from the loop

%xfer 
pi_clxfer = feedback(pi_design_K*pi_design_D*ems_xfer,sen_xfer*Kf*Df);  %feedback(G,H)

Q2.K = pi_design_K;
Q2.Z = pi_design_Z;
Q2.X = pi_clxfer;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%Q3: PID Dynamics

%PID controller parameters given in the question
pid_pole = -2*CF;
pid_zero1 = ems_p1;     %closest pole to the jw-axis
pid_zero2 = ems_p2;     %the other closest pole to the jw-axis

%PID controller gains
pid_Ki = 1;
pid_Kp = 1/pid_pole-(pid_zero1+pid_zero2)/(pid_zero1*pid_zero2);    %Kp = 1/p - (z1+z2)/(z1*z2)
pid_Kd = 1/(pid_zero1*pid_zero2)+pid_Kp/pid_pole;                   %Kd = (z-p)/(z*p), (finite difference derrivative), Fdd = -(p*s)/(s-p) [usually p = 2*CF]

%PID controller dynamics
pid_D = pid_Kp+pid_Ki*(1/s)+pid_Kd*(-pid_pole*s)/(s-pid_pole);      %D = [p(s-z)]/[z(s-p)] = Kp + Kd * Fdd + Ki * 1/s

Q3.Kp = pid_Kp;
Q3.Kd = pid_Kd;     
Q3.D = pid_D;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%Q4: Real Zeros

%find the ultimate gain (kappa)
pid_Ku = margin(Q3.D*ol_xfer);      %margin(DGH)

%xfer when K = Ku(kappa)/2
pid_K50 = pid_Ku*50/100;
pid_xfer50 = feedback(pid_K50*Q3.D*ems_xfer,sen_xfer*Df*Kf);    %feedback(G,H), closed loop

Q4.Ku = pid_Ku;
Q4.X = pid_xfer50;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%Q5: Complex Zeros

%PID controller parameters change to complex zeros from assignment 6
pid_zero1_comp = ems_p3;
pid_zero2_comp = ems_p4; 

%PID controller gains
pid_Ki_comp = 1;
pid_Kp_comp = 1/pid_pole-(pid_zero1_comp+pid_zero2_comp)/(pid_zero1_comp*pid_zero2_comp);   %Kp = 1/p - (z1+z2)/(z1*z2)
pid_Kd_comp = 1/(pid_zero1_comp*pid_zero2_comp)+pid_Kp_comp/pid_pole;                       %Kd = (z-p)/(z*p), (finite difference derrivative), Fdd = -(p*s)/(s-p) [usually p = 2*CF]

%PID controller dynamics
pid_D_comp = pid_Kp_comp+pid_Ki_comp*(1/s)+pid_Kd_comp*(-pid_pole*s)/(s-pid_pole);          %D = [p(s-z)]/[z(s-p)] = Kp + Kd * Fdd + Ki * 1/s

%find the ultimate gain
pid_Ku_comp = margin(pid_D_comp*ol_xfer);   %margin(DGH)

%xfer when K = Ku(kappa)/2
pid_K50_comp = pid_Ku_comp*50/100;
pid_xfer50_comp = feedback(pid_K50_comp*pid_D_comp*ems_xfer,sen_xfer*Df*Kf);    %feedback(G,H), closed loop

Q5.Ku = pid_Ku_comp;
Q5.X = pid_xfer50_comp;

%%%TUNING PARAMETERS
PID_CNTRL = 0.0690+0.630*(1/s)+0.00095*(2*CF)*s/(s+2*CF);   %adjust the numbers and check the stepinfo to meet the requirements
overall_xfer = feedback(1*PID_CNTRL*ems_xfer,sen_xfer*Df*Kf);
figure(20); clf;
step(overall_xfer);
stepinfo(overall_xfer)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%Q6: Tuned PID Control
%settle time: 0.1603
Q6.K = 0.630;           %K = Ki from PID_CNTRL because we assume Ki = 1
Q6.Z = [-46 -10.8];     %pzmap(PID_CNTRL), hard code the zeros obtained from D
Q6.X = overall_xfer;

%%%TEST PLOTS
if 0
    figure(10); clf;
    step(Q4.X);
    hold on; 
    step(Q5.X);
    grid on;
end
