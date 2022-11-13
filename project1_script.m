%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% MATLAB PROJECT 1 - ELEC 341 - IDIL BIL
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
Q2.Ga = minreal(Kdc*(WnTr_average^2)/(s^2+2*zeta*WnTr_average*s+WnTr_average^2));   %(V/V) xfer = Kdc*omega^2/(s^2+2*zeta*omega*s+omega^2)

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

%variable declaration
G1 = 100/(s+A*500);
G2 = 100/(s+B*600);
G3 = 10^5/(s+C*700);
G4 = 1500/(s+D*800);
H1 = 3/(s+F*5);

%transfer function found by simplifying the block diagram
tf_q3 = (G2*(G1+(G3/(1+G4*G3*H1))))/(1+((G2*G3*H1)/(1+G4*G3*H1)));

%gain is the final value when the input is 1/s so can be found by evaluating the transfer function at s = 0
temp_q3 = minreal(tf_q3);
[num den] = tfdata(temp_q3);
syms s;                             %define s back as a symbol
snum = poly2sym(cell2mat(num),s);   %polynomial to symbolic
sden = poly2sym(cell2mat(den),s);   %polynomial to symbolic
dcgain_num = subs(snum,s,0);        %s = 0 for the gain
dcgain_den = subs(sden,s,0);        %s = 0 for the gain
s = tf('s');                        %define s for laplace again
Q3.Kdc = double(dcgain_num/dcgain_den); 
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
Jr = F/15*10^(-7);      %(kg-m^2)
Br = G/30*10^(-6);      %(Nms)
Js = H/5*10^(-7);       %(kg-m^2)
Ms = A/4*10^(-3);       %(kg)
Bs = B/3;               %(Ns/m)
Ns = 3*10^(-2)/(2*pi);  %(m/rad) 

%variable declarations for mechanism and controller
Jf = C/3*10^(-7);        %(kg-m^2)
Bf = D/50;               %(Nms)
Nf = 10*10^(2)*(pi/180); %(rad/m)
Bt = E*10^(-3);          %(Nms)
Kt = F*30*10^(-3);       %(Nm)
L6 = 100*10^(-3);        %(m)
CF = 200;                %(Hz)

%nV = I/n * Z
Jp = Jr+Js+Mh*Ns^2+Ms*Ns^2+Jf*(Ns^2*Nf^2);  %simplification from the diagram 
Bp = Br+Bs*Ns^2+Bf*(Ns^2*Nf^2);             %simplification from the diagram
Ktp = Kt*(Ns^2*Nf^2);                       %calculated seperately cause connected to the load
Btp = Bt*(Ns^2*Nf^2);                       %calculated seperately cause connected to the load

%[sqr; swr; sIw] = A * [qr; wr; Iw] + B * [Vin; ](inputs)
A4 = [0 1 0;
     -Ktp/Jp (-Bp-Btp)/Jp Km/Jp;
     0 -Km/Lw -Rw/Lw];
B4 = [0 0 1/Lw].';

%[qr; ff] = C * [qr; wr; Iw] + B * [Vin; ](inputs)
C4 = [1 0 0;
     Ktp/L6 Btp/L6 0];
D4 = [0 0].';           %same as D = [0; 0]

ID = eye(3);
STM = inv(s*ID-A4);     %system transition matrix = (s*identity matrix - A)^(-1)
tf_q4 = C4*STM*B4+D4;   %transfer function = C*STM*B+D

Q4.Gj = minreal(tf_q4(1)*180/pi);    %with conversion from radians to degrees
Q4.Gt = minreal(tf_q4(2)/0.00981);   %with conversion from N to grams

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
Q5.Kfb = 1/(10^(-3)*Q3.Kdc);  %mV to volts

% looking at pzmap for sensor xfer (Q3.Hs) reveals that most dominant pole (rightmost x) is 75 rad/s
Nhat = 1.72;                        %CF*2*pi/[most dominant pole (rightmost x)*10] closest value to this on the N-Nhat table
Q5.Du = (CF*2*pi)/(Nhat*s+CF*2*pi); %(Hz to rad) CF/(Nhat+CF(Hz to rad))

looptf = Q2.Ga*Q4.Gj*Q3.Hs*10^(-3)*Q5.Du*Q5.Kfb; %open loop transfer function, calculated using the simulink module
Q5.GH = minreal(looptf);

%plots to check Q5
if 0
    figure(5); clf;
    step(Q5.GH);    %check is it gets stable in a reasonable amount of time
    grid on;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%Q6: Joint-Space Function

%in the simulink diagram
G6 = 1*Q2.Ga*Q4.Gj;
H6 = Q3.Hs*10^(-3)*Q5.Du*Q5.Kfb;
Q6.Xj  = minreal(G6/(1+G6*H6));    %G/(1+GH)

%plots to check Q6
if 0
    figure(6); clf;
    step(Q6.Xj);    %check is it gets stable in a reasonable amount of time
    grid on;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%Q7: Task-Space Function

%in the simulink diagram
G7 = 1*Q2.Ga;
H7 = Q4.Gj*Q3.Hs*10^(-3)*Q5.Du*Q5.Kfb;
Q7.Xt = minreal(Q4.Gt*G7/(1+G7*H7));    %G/(1+GH)

%plots to check Q7
if 0
    figure(7); clf;
    step(Q7.Xt);    %check is it gets stable in a reasonable amount of time
    grid on;
end
