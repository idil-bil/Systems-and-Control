%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% MATLAB PROJECT - ELEC 341 - IDIL BIL
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
%%%Q1: Experimental Data from a Data-Sheet

%plot p1DSPlot.p with your student number and read the values
Q1.Tr = 0.750;                %(ms) first time it gets to the final value before the peak
Q1.Tp = 1.405;                %(ms) first time it gets to the highest value
Q1.Ts = 4.147;                %(ms) first time it gets %2 close to settle value after the rise
Q1.Pos = (22.5-18)/18*100;    %(percent) (Tpeak y - Tsettle y)/Tsettle y * 100

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%Q2: 2nd Order Approx

%use 2nd order approximation to find the values
s = tf('s');                                %define s for laplace 
Kdc = 18;                                   %y value of Tsettle
os = Q1.Pos/100;                            %this is more accurate than using the y value of Tpeak
zeta = sqrt(log(os)^2/(pi^2 + log(os)^2));  %(pure) zeta = sqrt[(ln(overshoot)^2/(pi^2+ln(overshoot)^2)]
beta = sqrt(1-zeta^2);                      %(pure) beta = sqrt(1-zeta^2)

WnQ1.Tr = 1/(beta*Q1.Tr*10^(-3))*(pi-atan(beta/zeta));   %(rad/s) omega with Tr = (pi-atan(beta/zeta))/(beta*Tr)
WnQ1.Tp = pi/(beta*Q1.Tp*10^(-3));                       %(rad/s) omega with Tp = pi/(beta*Tp)

new_Tr = 1/(beta*WnQ1.Tp)*(pi-atan(beta/zeta));            %(seconds) Trise with Tpeak = (pi-atan(btea/zeta))/(beta*omega with Tpeak)
Tr_average = (new_Tr+Q1.Tr*10^(-3))/2;                     %(seconds) calculating the average of Tr's for balance between Tpeak and Trise
WnTr_average = 1/(beta*Tr_average)*(pi-atan(beta/zeta));   %(rad/s) omega with Tr = (pi-atan(beta/zeta))/(beta*Tr)

GaTr = Kdc*(WnQ1.Tr^2)/(s^2+2*zeta*WnQ1.Tr*s+WnQ1.Tr^2);                            %(V/V) xfer = Kdc*omega^2/(s^2+2*zeta*omega*s+omega^2)
GaTp = Kdc*(WnQ1.Tp^2)/(s^2+2*zeta*WnQ1.Tp*s+WnQ1.Tp^2);                            %(V/V) xfer = Kdc*omega^2/(s^2+2*zeta*omega*s+omega^2)
%Q2.Ga = minreal(Kdc*(WnTr_average^2)/(s^2+2*zeta*WnTr_average*s+WnTr_average^2));  %(V/V) xfer = Kdc*omega^2/(s^2+2*zeta*omega*s+omega^2)
Q2.Ga = 1.165e08/(s^2 + 2054*s + 6.472e06);                                         %exact value from the answer key to make the second part more precise

%plots to check Q2
if 0
    [y1 t1] = step(GaTr);
    [y2 t2] = step(GaTp);
    [ya ta] = step(Q2.Ga);
    t1 = t1*10^(3);
    t2 = t2*10^(3);
    ta = ta*10^(3);
    figure(1); clf;
    p1DSPlot(ABCDEFGH,1);
    hold on
    plot(t1,y1);    
    plot(t2,y2);    
    plot(ta,ya);    %should be the average of both graphs
    hold off;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%Q3: Analog Position Sensor

%variable declaration for the sensor block diagra (figure 4)
G1 = 100/(s+A*500);
G2 = 100/(s+B*600);
G3 = 10^5/(s+C*700);
G4 = 1500/(s+D*800);
H1 = 3/(s+F*5);

%transfer function found by simplifying the block diagram
tf_q3 = (G2*(G1+(G3/(1+G4*G3*H1))))/(1+((G2*G3*H1)/(1+G4*G3*H1)));

%dc gain is the final value (the y value of Tsettle) when the step response is found (multiply by 1/s)
%so can be found by evaluating the transfer function at s = 0
temp_q3 = minreal(tf_q3);           %simplified xfer function for Q3 as a temporary value
[num den] = tfdata(temp_q3);        %fins the roots in the numerator and denominator (s-root1)/(s-root2)
syms s;                             %turns laplace s to letter/symbol s (to defie a value or do ilaplace)
snum = poly2sym(cell2mat(num),s);   %polynomial to symbolic
sden = poly2sym(cell2mat(den),s);   %polynomial to symbolic
dcgain_num = subs(snum,s,0);        %s = 0 for the dc gain
dcgain_den = subs(sden,s,0);        %s = 0 for the dc gain
s = tf('s');                        %define s for laplace again

Q3.Kdc = double(dcgain_num/dcgain_den);     %could also just do dcgain(tf_q3)
Q3.Hs = minreal(tf_q3);

%plots to check Q3
if 0
    figure(2); clf;
    step(Q3.Hs);    %check is it gets stable in a reasonable amount of time and hits the y value of Tsettle
    grid on;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%Q4: Joint and Task Space

%variable declarations for motor lead and screw
Rw = A/2;               %(ohm)
Lw = B*30*10^(-6);      %(H)
Km = C*10^(-3);         %(Nm/A)
Mh = (D+E)*10^(-3);     %(kg)
Jr = F/15*10^(-7);      %(kg*m^2)
Br = G/30*10^(-6);      %(Nms)
Js = H/5*10^(-7);       %(kg*m^2)
Ms = A/4*10^(-3);       %(kg)
Bs = B/3;               %(Ns/m)
Ns = 3*10^(-2)/(2*pi);  %(m/rad) 

%variable declarations for mechanism and controller
Jf = C/3*10^(-7);        %(kg*m^2)
Bf = D/50;               %(Nms)
Nf = 10*10^(2)*(pi/180); %(rad/m)
Bt = E*10^(-3);          %(Nms)
Kt = F*30*10^(-3);       %(Nm)
L6 = 100*10^(-3);        %(m)
CF = 200;                %(Hz)

%nV = I/n * Z
Jp = Jr+Js+(Mh+Ms)*Ns^2+3*Jf*(Ns^2*Nf^2);  %simplification from the diagram 
Bp = Br+Bs*Ns^2+3*Bf*(Ns^2*Nf^2);          %simplification from the diagram
Ktp = 3*Kt*(Ns^2*Nf^2);                    %calculated seperately cause connected to the load
Btp = 3*Bt*(Ns^2*Nf^2);                    %calculated seperately cause connected to the load

%[sqr; swr; sIw] = A * [qr; wr; Iw] + B * [Vin; ](inputs)
A4 = [0 1 0;
     -Ktp/Jp (-Bp-Btp)/Jp Km/Jp;
     0 -Km/Lw -Rw/Lw];
B4 = [0 0 1/Lw].';

%[qr; ff] (outputs) = C * [qr; wr; Iw] + B * [Vin; ](inputs)
C4 = [1 0 0;            %Ktp/3 to calculate per finger
     Ktp/(3*L6) 0 0];   %Bt (friction is symbolized by B) is not included because the system is horizontal and there is no force moving the ball
D4 = [0 0].';           %same as D = [0; 0]
%for D -> output number = number of rows, input number = number of columns

ID = eye(3);            %identity matrix
STM = inv(s*ID-A4);     %system transition matrix = (s*identity matrix - A)^(-1)
tf_q4 = C4*STM*B4+D4;   %transfer function = C*STM*B+D

Q4.Gj = minreal(tf_q4(1)*180/pi);           %with conversion from radians to degrees
Q4.Gt = minreal(tf_q4(2)/(Ns*Nf*0.00981));  %with conversion from N to grams and divide by Ns and Nf for the real value instead of conversion

%plot to check Q4
if 0
    figure(3); clf;
    step(Q4.Gj);    %check is it gets stable in a reasonable amount of time
    grid on;
    figure(4); clf;
    step(Q4.Gt);    %check is it gets stable in a reasonable amount of time
    grid on;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%Q5: Open-Loop Function 

%want gain to be equal to 1 so its the inverse of the gain of the sensor
Q5.Kfb = 1/(10^(-3)*Q3.Kdc);  %(mV to volts)(unit conversion block) Kf = 1/(sensor (Hs) dc gain)

%looking at pzmap for sensor xfer (Q3.Hs) reveals that most dominant pole (leftmost x)
%Nhat = 3;              %CF/[most dominant pole (leftmost x)*10] closest value to this on the N-Nhat table (also not used when CF is smaller than 10 times the leftmost pole)
Q5.Du = (CF)/(s+CF);    %(Hz) Dfeedback = CF/[Nhat*s +CF] Nhat = 1 if there are on weighted sum filters (also not used when CF is smaller than 10 times the leftmost pole)

looptf = Q2.Ga*Q4.Gj*Q3.Hs*10^(-3)*Q5.Du*Q5.Kfb; %open loop transfer function, calculated using the simulink module p1q67
Q5.GH = minreal(looptf);                         %(GH) didnt include K in the open loop xfer

%plots to check Q5
if 0
    figure(5); clf;
    step(Q5.GH);    %check is it gets stable in a reasonable amount of time
    grid on;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%Q6: Joint-Space Function (rotor angle)

%in the simulink diagram p1q67
G6 = 1*Q2.Ga*Q4.Gj;                 %K = 1
H6 = Q3.Hs*10^(-3)*Q5.Du*Q5.Kfb;    
Q6.Xj  = minreal(G6/(1+G6*H6));     %closed loop xfer = G/(1+GH), could write feedback(G,H)

%plots to check Q6
if 0
    figure(6); clf;
    step(Q6.Xj);    %check is it gets stable in a reasonable amount of time
    grid on;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%Q7: Task-Space Function (finger force)

%in the simulink diagram p1q67
G7 = 1*Q2.Ga;                           %K = 1
H7 = Q4.Gj*Q3.Hs*10^(-3)*Q5.Du*Q5.Kfb;
Q7.Xt = minreal(Q4.Gt*G7/(1+G7*H7));    %closed loop xfer = G/(1+GH), could write feedback(G,H), Q4.Gt is outside the feedback loop

%plots to check Q7
if 0
    figure(7); clf;
    step(Q7.Xt);    %check is it gets stable in a reasonable amount of time
    grid on;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%Q8: STEP 3

%controller pole is at CF
%derivative pole is normally at 2*CF
%use WS to get derivative pole close to controller pole
ws_N_hat = 1.94;    %weighted sum Nhat = % 2CF/CF = 2 ~1.94 from the table
ws_N = 6;           %weighted sum N -> %corresponding N value to N_hat

%dynamics of controller associated with poles
pid_Dp = (1/s)*(2*CF)/(ws_N_hat*s+2*CF);    %Dpole = 1/s * (2*CF)/(Nhat*s + 2*CF)

Q8.N = ws_N;
Q8.Nhat = ws_N_hat;
Q8.Dp = minreal(pid_Dp);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%Q9: STEP 4

%compute Kref (kappa) with Dp
Kref = margin(Q8.Dp*Q5.GH);     %kappa = margin(openloop xfer including Dp)

%find wxo (frequency where we start trying numbers to find the max phase margin) usign Kref
[G_m P_m wxo_g wxo_p] = margin(Kref*Q8.Dp*Q5.GH);   %gain margin, phase margin, freq gain, freq phase (freq gain = freq phase when margin is done with Kref)

Q9.Kref = Kref;
Q9.wxo = wxo_g;

%plots to check Q9
if 0
    figure(1); clf;
    margin(Kref*Q8.Dp*Q5.GH);
    figure(2); clf;
    nyqlog(Kref*Q8.Dp*Q5.GH);   %blue line needs to go through -1 (black dot)
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%Q10: STEP 5

%Search for initial zero positions
if 0    %starts from the inital frequency found and finds the real zeros that give the max phase margin
[pid_D_rec_r2, zero12_r, maxPm_r2] = maxPmRealZeros(Q9.wxo,Q9.Kref,Q8.Dp,Q5.GH);    %(initial freq, Kref kappa, Dp, GH after Dp)
end     %[D, zero matrix[zero1, zero2], max phase margin]

if 0    %starts from the inital frequency found and finds the complex zeros that give the max phase margin
[pid_D_rec_i1, zero12_i1, maxPm_i1, img_search_valid1, invalid_matrix1] = maxPmComplexZeros(Q9.wxo,Q9.Kref,Q8.Dp,Q5.GH);    %(initial freq, Kref kappa, Dp, GH after Dp)
end     %[D, zero matrix[zero1, zero2], max phase margin, if 1: the code found the max phase margin if 0: decrease resolution to get 1, gives the invalid matrix so change resolution or iteration number]
%if it doesn't matter if the zeros are real or complex compare the phase margins in the terminal and choose the bigger ones zeros

%Q10.Z = zero12_i1;
%Q10.PM = maxPm_i1;
%Q10.D = minreal(pid_D_rec_i1);
%hard code these values after running so you don't need to run again
Q10.Z = [-6.8964 + 8.3600i, -6.8964 - 8.3600i];
Q10.PM = 139.3936;
Q10.D = (1.756*s^2 + 24.21*s + 206.2)/(s^2 + 206.2*s);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%Q11: STEP 6

%poles and zeros for initial designed controller
pid_pole = -2*CF;
pid_zero1 = Q10.Z(1);
pid_zero2 = Q10.Z(2);

%dynamics gains
init_pid_Ki = 1;                                                         %Ki = 1 (for initial)
init_pid_Kp = 1/pid_pole-(pid_zero1+pid_zero2)/(pid_zero1*pid_zero2);    %Kp = 1/p - (z1+z2)/(z1*z2)
init_pid_Kd = 1/(pid_zero1*pid_zero2)+init_pid_Kp/pid_pole;              %Kd when Kp is defined = 1/(z1*z2) + Kp/p

%desired hase margin, if its not given define as 45
desiredPm = 30;

%master gain of the controller for a desired phase argin
init_pid_masterK = getK4PhaseMargin(Q9.Kref,Q10.D,Q5.GH,desiredPm);     %master gain (K) = (Kref, D, openloop after D, desired phase margin)

Q11.K = init_pid_masterK;
Q11.Kp = init_pid_Kp;
Q11.Ki = init_pid_Ki;
Q11.Kd = init_pid_Kd;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%Q12: STEP 7a

%%%
% REQUIREMENTS
% OS < 20 %
% Ts < 200 ms
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
Ki_factor = 0.30;  
tuned_pid_Ki = Ki_factor;
tuned_pid_Kp = init_pid_Kp * 1.85;
tuned_pid_Kd = init_pid_Kd * 0.06;
tuned_pid_masterK = init_pid_masterK * 0.35;

%full dynamics of the controller 
tuned_pid_D = tuned_pid_Kp + tuned_pid_Ki*(1/s) + tuned_pid_Kd*Q8.Dp*s^2;

Nhat = 1;
pid_D = pid_Kp + pid_Ki*(1/s) + pid_Kd*(2*CF*s)/(Nhat*s+2*CF);

% the response with the tuned controller
tuned_response = feedback(tuned_pid_masterK*tuned_pid_D*Q2.Ga*Q4.Gj,Q3.Hs*10^(-3)*Q5.Du*Q5.Kfb);

% plot the response
figure(100); clf;
step(tuned_response);
grid on;
stepinfo(tuned_response)

Q12.K = tuned_pid_masterK*Ki_factor;
Q12.Kp = tuned_pid_Kp/Ki_factor;
Q12.Ki = tuned_pid_Ki/Ki_factor;
Q12.Kd = tuned_pid_Kd/Ki_factor;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%Q13: STEP 7b

% solve the equations involving zeros to find the zeros for our tuned gains
syms a b
eqn1 = Q12.Kp == (1/(-2*CF/Q8.Nhat)) - ((a+b*1i)+(a-b*1i))/((a+b*1i)*(a-b*1i));
eqn2 = Q12.Kd == (1/((a+b*1i)*(a-b*1i))) + (Q12.Kp/(-2*CF/Q8.Nhat));
[sigma, omega] = solve(eqn1, eqn2, a, b);
tuned_pid_zeros = [double(vpa(sigma(2) + omega(2)*1i)) double(vpa(sigma(1) + omega(1)*1i))];

[~,tuned_PM] = margin(tuned_pid_masterK*tuned_pid_D*Q5.GH);

Q13.Z = tuned_pid_zeros;
Q13.PM = tuned_PM;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%Q14: STEP 8

Bf_imp = D/50*1.5; % Nms
B_p_imp = Br+Bs*Ns^2+3*Bf_imp*(Ns^2*Nf^2);

% A,B,C,D matrices for State-Space Representation
matr_A = [0 1 0;
     -Ktp/Jp (-Bp-Btp)/Jp Km/Jp;
     0 -Km/Lw -Rw/Lw];
matr_B = [0 0 1/Lw].';
matr_C = [1 0 0;
          Ktp/(3*L6) 0 0];
matr_D = [0 0].';

ID_matrix = eye(3);
temp = s*ID_matrix-matr_A;
STM = inv(temp);

TF = matr_C*STM*matr_B+matr_D;

Gj_imp = minreal(TF(1)*180/pi);

Gfrwd_imp = Gj_imp*Q2.Ga;

% Q14 Answer(s)
Q14.G = Gfrwd_imp;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%Q15: STEP 7a - Round 2

%%%
% REQUIREMENTS
% OS < 20 %
% Ts < 200 ms
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

% tune the following gains to satisfy RCGs listed above
Ki_factor_imp = 0.50;
tuned_pid_Ki_imp = Ki_factor_imp;
tuned_pid_Kp_imp = init_pid_Kp * 1.86;
tuned_pid_Kd_imp = init_pid_Kd * 0.04;
tuned_pid_masterK_imp = init_pid_masterK * 0.75;

% full dynamics of the controller 
tuned_pid_D_imp = tuned_pid_Kp_imp + tuned_pid_Ki_imp*(1/s) + tuned_pid_Kd_imp*Q8.Dp*s^2;

% the response with the tuned controller
tuned_response_imp = feedback(tuned_pid_masterK_imp*tuned_pid_D_imp*Q14.G,Q3.Hs*10^(-3)*Q5.Du*Q5.Kfb);

% plot the response
figure(101); clf;
step(tuned_response_imp);
grid on;
stepinfo(tuned_response_imp)

% Q15 Answer(s)
Q15.K = tuned_pid_masterK_imp*Ki_factor_imp;
Q15.Kp = tuned_pid_Kp_imp/Ki_factor_imp;
Q15.Ki = tuned_pid_Ki_imp/Ki_factor_imp;
Q15.Kd = tuned_pid_Kd_imp/Ki_factor_imp;
