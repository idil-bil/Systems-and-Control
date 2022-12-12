%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% MATLAB ASSIGNMENT 7 - ELEC 341 - IDIL BIL
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
rate_ol_xfer = ems_xfer*sen_xfer*Df*Kf;     %K is not inluded since its not calculated yet  (GH)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%Q1: Lead Dynamics

%cancel out the rightmost by putting a zero on top of the pole (most dominant, closest to jw-axis) (ems_p1)
%given that FDD has pole at -2CF
speed_lead_dyn = (s-ems_p1)/(s+2*CF);           %(s-z)/(s-p)
speed_lead_gain = 2*CF/(abs(ems_p1));           %p/z
speed_lead_KD = (ems_p1+2*CF)/abs(ems_p1*2*CF); %Kd = (z-p)/(z*p), (finite difference derrivative) Fdd = -(p*s)/(s-p) [usually p = 2*CF]

Q1.Kd = speed_lead_KD;
Q1.D = speed_lead_gain*speed_lead_dyn;          %D = [p(s-z)]/[z(s-p)] = Kp + Kd * Fdd + Ki * 1/s

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%Q2: Steady State Error %

%find the ultimate gain Kappa
speed_lead_kappa = margin(Q1.D*rate_ol_xfer);   %DGH

%steady state error: K = Ku(kappa, ultimate gain) * 50%
speed_lead_k50 = speed_lead_kappa/2;
speed_lead_xfer50 = feedback(speed_lead_k50*Q1.D*ems_xfer,sen_xfer*Kf*Df);  %feedback(G,H)
speed_lead_Ess50 = abs(1-dcgain(speed_lead_xfer50))*100;                    %Ess = abs(1-dcgain(Closed loop function))*100

%steady state error: K = Ku(kappa, ultimate gain) * 99%
speed_lead_k99 = speed_lead_kappa*99/100;
speed_lead_xfer99 = feedback(speed_lead_k99*Q1.D*ems_xfer,sen_xfer*Kf*Df);  %feedback(G,H)
speed_lead_Ess99 = abs(1-dcgain(speed_lead_xfer99))*100;                    %Ess = abs(1-dcgain(Closed loop function))*100

Q2.Ku = speed_lead_kappa;
Q2.Ess = speed_lead_Ess50;
Q2.Ess99 = speed_lead_Ess99;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%Q3: Tuned Rate Control

if 0                            %for minimum ess
    minimum_Ess = 100;

    for cnt1 = 0:-1:-100        %from 0 to -100 by -1
        for cnt2 = 1:1:95
            valid = 1;
            testZ = cnt1;       %will get incremented

            %adjust new xfer for test
            speed_lead_dyn_test = (s-testZ)/(s+2*CF);
            speed_lead_gain_test = 2*CF/abs(testZ); 
            speed_lead_D_test = speed_lead_gain_test*speed_lead_dyn_test;

            %find new kappa value
            speed_lead_kappa_test = margin(speed_lead_D_test*rate_ol_xfer);

            %try K values
            speed_lead_K_test = speed_lead_kappa_test*cnt2/100;

            %get CL xfer
            speed_lead_clxfer_test = feedback(speed_lead_K_test*speed_lead_D_test*ems_xfer,sen_xfer*Kf*Df);

            %check peak value
            [ytest xtest] = step(speed_lead_clxfer_test);

            for cnt3 = 1:length(ytest)
                if(ytest(cnt3) >= 1.2)      %overshoot requirement is 20%
                    valid = 0;
                    break;
                end
            end

            %check min Ess if the ess requirement is not met
            speed_lead_Ess_test = abs(1-dcgain(speed_lead_clxfer_test))*100;

            if(valid == 1)
                if((speed_lead_Ess_test < minimum_Ess) && (speed_lead_Ess_test < Q2.Ess))   %&& is to compare to the ess of Q2 as the quesition asked
                    minimum_Ess = speed_lead_Ess_test;
                    speed_design_K = speed_lead_K_test;
                    speed_design_K_percent = cnt2;
                    speed_design_Z = testZ;
                end
            end
        end
    end
end

%hard coded design values
speed_design_K = 0.1965;    %get this values from the loop
speed_design_Z = -39;       %get this values from the loop

speed_design_D = 2*CF/abs(speed_design_Z)*(s-speed_design_Z)/(s+2*CF);      %D = [p(s-z)]/[z(s-p)] = Kp + Kd * Fdd + Ki * 1/s

Q3.K = speed_design_K;
Q3.Z = speed_design_Z; 
Q3.X = feedback(speed_design_K*speed_design_D*ems_xfer,sen_xfer*Df*Kf);     %feedback(G,H)

speed_design_Ess = abs(1-dcgain(Q3.X))*100;                                 %Ess = abs(1-dcgain(Closed loop function))*100

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%Q4: Ultimate Gain

% OL xfer from Prev. Assignment
pos_ol_xfer = ems_xfer*(1/s)*sen_xfer*Df*Kf;

% find ultimate gain
pos_lead_kappa = margin(Q1.D*pos_ol_xfer);

% K = Ku/2
pos_lead_k50 = pos_lead_kappa/2;

% find the CL xfer
pos_lead_xfer50 = feedback(pos_lead_k50*Q1.D*ems_xfer*(1/s),sen_xfer*Kf*Df);

Q4.Ku = pos_lead_kappa;
Q4.X = pos_lead_xfer50;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%Q5: Tuned Position Control
minimum_Ts = 10000000;

if 0                        %for minimum settle time
    for cnt1 = 0:-1:-100
        for cnt2 = 1:1:95
            valid = 1;
            testZ = cnt1;   %will get incremented

            %adjust new xfer for test
            pos_lead_dyn_test = (s-testZ)/(s+2*CF);
            pos_lead_gain_test = 2*CF/abs(testZ); 
            pos_lead_D_test = pos_lead_gain_test*pos_lead_dyn_test;

            %find new kappa value
            pos_lead_kappa_test = margin(pos_lead_D_test*pos_ol_xfer);

            %try K values
            pos_lead_K_test = pos_lead_kappa_test*cnt2/100;

            %get CL xfer
            pos_lead_clxfer_test = feedback(pos_lead_K_test*pos_lead_D_test*ems_xfer*(1/s),sen_xfer*Kf*Df);

            %check peak value
            [ytest xtest] = step(pos_lead_clxfer_test);

            for cnt3 = 1:length(ytest)
                if(ytest(cnt3) >= 1.2)      %overshoot requirement is 20%
                    valid = 0;
                    break;
                end
            end

            % if(stepinfo(pos_lead_clxfer_test).SettlingTime < 0.2)     %if the requirement was that the settling time was smaller than 200 ms
            %     valid = 0;
            % end

            %check min Ts

            if(valid == 1)
                if(stepinfo(pos_lead_clxfer_test).SettlingTime < minimum_Ts)        %change the "SettlingTime" to "Overshoot" if you want to find min overshoot
                    minimum_Ts = stepinfo(pos_lead_clxfer_test).SettlingTime;
                    pos_design_K = pos_lead_K_test;
                    pos_design_K_percent = cnt2;
                    pos_design_Z = testZ;
                end
            end
        end
    end
end

%hard coded design values
pos_design_K = 0.4412;  %get this value from the loop
pos_design_Z = -10;     %get this value from the loop

pos_design_D = 2*CF/abs(pos_design_Z)*(s-pos_design_Z)/(s+2*CF);              %D = [p(s-z)]/[z(s-p)] = Kp + Kd * Fdd + Ki * 1/s

Q5.K = pos_design_K;
Q5.Z = pos_design_Z;
Q5.X = feedback(pos_design_K*pos_design_D*ems_xfer*(1/s),sen_xfer*Kf*Df);     %feedback(G,H)

%%%TEST PLOTS
if 0
    figure(1); clf;
    pzmap(rate_ol_xfer);
    figure(2); clf;
    pzmap(Q1.D*ems_xfer*sen_xfer*Kf*Df);
end

if 0
    figure(3); clf;
    step(speed_lead_xfer50);
    figure(4); clf;
    step(speed_lead_xfer99);
end

if 0
    figure(5); clf;
    step(Q3.X);
end

if 0
    figure(6); clf;
    step(Q5.X);
end