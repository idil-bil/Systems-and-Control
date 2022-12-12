%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% MATLAB ASSIGNMENT 6 - ELEC 341 - IDIL BIL
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%Q1: Rate Control

%pole-zero declarations for EMS (G - controller) and SEN (H - feedback) (figure 2)
ems_p1 = -A;            %X
ems_p2 = -3*B;          %X
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
[rate_Gm rate_Pm] = margin(rate_ol_xfer);   %[GainMargin PhaseMargin] = margin(open loop xfer except controller (GH))
rate_Kappa = rate_Gm;                       %ultimate (maximum) gain Ku (kappa) = Gain margin
rate_K = rate_Kappa/2;                      %P controller K (give as half of Kappa in the question)

%closed loop xfer
rate_cl_xfer = feedback(rate_K*ems_xfer,sen_xfer*Df*Kf);    %feedback equation = G/(1+GH)

%calculations for steady state error (used an easier way in the later assignments)
rate_KGH = rate_K*rate_ol_xfer;                                         %actual open loop xfer since K is included
[num den] = tfdata(rate_KGH);
syms s;                                                                 %turns laplace s to letter s (to define a value or do ilaplace)
rate_KGH_sym_s = poly2sym(cell2mat(num),s)/poly2sym(cell2mat(den),s);
s = 0;                                                                  %s = 0 for the steady state error
rate_Ess = eval(1/(1+rate_KGH_sym_s))*100;                              %(%)steady state error = 1/(1+KGH)
s = tf("s");                                                            %define s for laplace

Q1.GH = rate_ol_xfer;
Q1.Ku = rate_Kappa;
Q1.X = rate_cl_xfer;
Q1.Ess = rate_Ess;

%plots to check Q1
if 0
    figure(1); clf;
    pzmap(Q1.GH);   %check if the poles and zeros are at the correct place
    grid on;        %make sure all the poles and zeros are on the left hand side
    figure(2); clf;
    impulse(Q1.GH);
    grid on;
    figure(3); clf;
    step(Q1.GH);    
    grid on;

    figure(4); clf;
    impulse(Q1.X);
    grid on;
    figure(5); clf;
    step(Q1.X);     %settle value should be 1-steady state error
    grid on;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%Q2: Speed Control

%open loop xfer
pos_ol_xfer = ems_xfer*(1/s)*sen_xfer*Df*Kf;    %K is not inluded since its not calculated yet
[pos_Gm pos_Pm] = margin(pos_ol_xfer);          %[GainMargin PhaseMargin] = margin(open loop xfer except controller)
pos_Kappa = pos_Gm;                             %ultimate (maximum) gain Ku (kappa) = Gain margin
pos_K = pos_Kappa/2;                            %P controller K (give as half of Kappa in the question)

%closed loop transfer
pos_cl_xfer = feedback(ems_xfer*pos_K*(1/s),sen_xfer*Df*Kf);    %feedback equation = G/(1+GH) 

%calculations for steady state error (used an easier way in the later assignments)
pos_KGH = pos_K*pos_ol_xfer;                                            %actual open loop xfer since K is included
[num den] = tfdata(pos_KGH);                                            
syms s;                                                                 %turns laplace s to letter s (to define a value or do ilaplace)
pos_KGH_sym_s = poly2sym(cell2mat(num),s)/poly2sym(cell2mat(den),s);
s = 0;                                                                  %s = 0 for the steady state error
pos_Ess = eval(1/(1+pos_KGH_sym_s))*100;                                %(%)steady state error = 1/(1+KGH)
s = tf("s");                                                            %define s for laplace

Q2.GH = pos_ol_xfer;
Q2.Ku = pos_Kappa;
Q2.X = pos_cl_xfer;
Q2.Ess = pos_Ess;

%plots to check Q2
if 0
    figure(6); clf;
    pzmap(Q2.GH);   %check if the poles and zeros are at the correct place
    grid on;        %make sure all the poles and zeros are on the left hand side
    figure(7); clf;
    impulse(Q2.GH);
    grid on;
    figure(8); clf;
    step(Q2.GH);
    grid on;

    figure(9); clf;
    impulse(Q2.X);
    figure(10); clf;
    grid on;
    step(Q2.X);     %settle value should be 1-steady state error
    grid on;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Q3: Rate RCGs

%Requirement (must), Constraint (limit), Goal (optimization) 1
rate_rcg1_ess = 0.3;                                                    %requirement = 30% steady state error
[num den] = tfdata(rate_ol_xfer);           
syms s;                                                                 %turns laplace s to letter s (to define a value or do ilaplace)
rate_GH_sym_s = poly2sym(cell2mat(num),s)/poly2sym(cell2mat(den),s);
s = 0;                                                                  %s = 0 for the steady state error
rate_GH = eval(rate_GH_sym_s);
rate_rcg1_K = (1/rate_rcg1_ess-1)/rate_GH;                              %1/(1+KGH) should give the RCG, rearranging gives (1/Ess-1)/GH = K

%Requirement (must), Constraint (limit), Goal (optimization) 2
%requirement = 0% overshoot, min steady state error (close to 1 when the step is plotted)
rate_rcg2_K1 = Q1.Ku/100*10; 
rate_rcg2_K2 = Q1.Ku/100*7; 
rate_rcg2_K3 = Q1.Ku/100*5; 
rate_rcg2_K4 = Q1.Ku/100*3; 
rate_rcg2_K5 = Q1.Ku/100*1; 
%change the percentage and check the step plots to see which one fits the requirements the best
if 0
    figure(10); clf;
    test = feedback(ems_xfer*rate_rcg2_K1,sen_xfer*Df*Kf); 
    step(test);
    hold on;
    test = feedback(ems_xfer*rate_rcg2_K2,sen_xfer*Df*Kf); 
    step(test);
    test = feedback(ems_xfer*rate_rcg2_K3,sen_xfer*Df*Kf); 
    step(test);
    test = feedback(ems_xfer*rate_rcg2_K4,sen_xfer*Df*Kf); 
    step(test);
    test = feedback(ems_xfer*rate_rcg2_K5,sen_xfer*Df*Kf); 
    step(test);
    grid on
end

Q3.K1 = rate_rcg1_K;
Q3.K2 = rate_rcg2_K3;   %really small number so 0 works too

%plots to check Q3
if 0
    figure(11); clf;
    rate_cl_xfer_rcg1 = feedback(ems_xfer*Q3.K1,sen_xfer*Df*Kf); 
    step(rate_cl_xfer_rcg1);
    grid on;
    figure(12); clf;
    rate_cl_xfer_rcg2 = feedback(ems_xfer*Q3.K2,sen_xfer*Df*Kf); 
    step(rate_cl_xfer_rcg2,0.5);
    grid on;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Q4: Position RCGs

%Requirement (must), Constraint (limit), Goal (optimization) 1
%requirement = 10% overshoot
pos_rcg1_K = Q2.Ku/100*?;   %change the percentage and check the step plots to see which one fits the requirements the best

%Requirement (must), Constraint (limit), Goal (optimization) 2
%requirement = 0% overshoot, min settle time
pos_rcg2_K = Q2.Ku/100*?;   %change the percentage and check the step plots to see which one fits the requirements the best

%plots to check Q4
if 0
    s = tf("s");
    figure(13); clf;
    pos_cl_xfer_rcg1 = feedback(ems_xfer*Q4.K1*(1/s),sen_xfer*Df*Kf);
    step(pos_cl_xfer_rcg1,2);
    figure(14); clf;
    pos_cl_xfer_rcg2 = feedback(ems_xfer*Q4.K2*(1/s),sen_xfer*Df*Kf);
    step(pos_cl_xfer_rcg2);
    grid on;
end

Q4.K1 = pos_rcg1_K;
Q4.K2 = pos_rcg2_K;     %really small number so 0 works too
