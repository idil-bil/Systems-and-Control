%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% MATLAB ASSIGNMENT 2 - ELEC 341 - IDIL BIL
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
%%%Q1: Experimental Data from the Data Sheet

%plot a2DSPlot(ABCDEFGH)
Q1.Tr = 0.140;  %(ms) first time it gets to the final value before the peak
Q1.Tp = 0.380;  %(ms) first time it gets to the highest value
Q1.Ts = 0.758;  %(ms) first time it gets %2 close to settle value after the rise
Q1.Pos = 42.75; %(percent) (Tpeak y - Tsettle y)/Tsettle y * 100

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%Q2: Approximate #1

os = Q1.Pos/100;                                     %this is more accurate than using the y value of Tpeak
Q2.Z = sqrt(log(os)^2/(pi^2 + log(os)^2));           %(pure) zeta = sqrt[(ln(overshoot)^2/(pi^2+ln(overshoot)^2)]
beta = sqrt(1-Q2.Z^2);                               %(pure) beta = sqrt(1-zeta^2)
Q2.Wn = 1/(beta*Q1.Tr*10^(-3))*(pi-atan(beta/Q2.Z)); %(rad/s) omega with Tr = (pi-atan(beta/zeta))/(beta*Tr)
Kdc = 24;                                            %gain = y value of Tsettle
s = tf('s');                                         %define s for laplace
Q2.G = Kdc*(Q2.Wn^2)/(s^2+2*Q2.Z*Q2.Wn*s+Q2.Wn^2);   %(V/V) xfer = Kdc*omega^2/(s^2+2*zeta*omega*s+omega^2)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%Q3: Approximate #2

Q3.Wn = pi/(beta*Q1.Tp*10^(-3));                    %(rad/s) omega with Tp = pi/(beta*Tp)
Q3.G = Kdc*(Q3.Wn^2)/(s^2+2*Q2.Z*Q3.Wn*s+Q3.Wn^2);  %(V/V) xfer = Kdc*omega^2/(s^2+2*zeta*omega*s+omega^2)

%plots to check Q3
if 0
    [y t] = step(Q3.G);
    t = t*10^3;
    figure(1); clf;
    a2DSPlot(ABCDEFGH)
    hold on
    plot(t,y);  %want our approximation to be as close to the original plot as possible
    hold off
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%Q4: Approximate #3
os4 = os*2/3;                                      %scale the overshoot as asked
Q4.Z = sqrt(log(os4)^2/(pi^2 + log(os4)^2));       %(pure) zeta = sqrt[(ln(overshoot)^2/(pi^2+ln(overshoot)^2)]
beta4 = sqrt(1-Q4.Z^2);                            %(pure) beta = sqrt(1-zeta^2)
Q4.Wn = pi/(beta4*Q1.Tp*10^(-3));                  %(rad/s) omega with Tp = pi/(beta*Tp) 
Q4.G = Kdc*(Q4.Wn^2)/(s^2+2*Q4.Z*Q4.Wn*s+Q4.Wn^2); %(V/V) xfer = Kdc*omega^2/(s^2+2*zeta*omega*s+omega^2)

%plots to check Q4
if 0
    [y t] = step(Q4.G);
    t = t*10^3;
    figure(1); clf;
    a2DSPlot(ABCDEFGH)
    hold on
    plot(t,y);  %want our approximation to be as close to the original plot as possible
    hold off
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%Q5: Approximate #4
os5 = os*1/3;                                       %scale the overshoot as asked
Q5.Z = sqrt(log(os5)^2/(pi^2 + log(os5)^2));        %(pure) zeta = sqrt[(ln(overshoot)^2/(pi^2+ln(overshoot)^2)]
beta5 = sqrt(1-Q5.Z^2);                             %(pure) beta = sqrt(1-zeta^2)
Q5.Wn = pi/(beta5*Q1.Tp*10^(-3));                   %(rad/s) omega with Tp = pi/(beta*Tp) 
Q5.G = Kdc*(Q5.Wn^2)/(s^2+2*Q5.Z*Q5.Wn*s+Q5.Wn^2);  %(V/V) xfer = Kdc*omega^2/(s^2+2*zeta*omega*s+omega^2)

%plots to check Q5
if 0
    [y t] = step(Q5.G);
    t = t*10^3;
    figure(1); clf;
    a2DSPlot(ABCDEFGH)  %want our approximation to be as close to the original plot as possible
    hold on
    plot(t,y);
    hold off
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%Q6 Approximate #5:

Q6.Z = 1;                                            %(pure) zeta = 1 when critically damped
beta6 = 0;                                           %(pure) beta = 0 when critically damped
%Q6.Wn = 4/(Q6.Z*Q1.Ts*10^(-3));                     %(rad/s) when critically damped Tr and Tp can't be used so omega = 4/(zeta*Ts)
%Q6.G = Kdc*(Q6.Wn^2)/(s^2+2*Q6.Z*Q6.Wn*s+Q6.Wn^2);  %(V/V) xfer = Kdc*omega^2/(s^2+2*zeta*omega*s+omega^2)

%More accurate way
syms s t w                                                                        %turns laplace s to letter s (to defie a value or do ilaplace)
step_response_s = Kdc/s*w^2/(s^2+2*w*s+w^2);
step_response_t = ilaplace(step_response_s,t);
Q6.Wn = double(vpasolve(subs(step_response_t,t,Q1.Ts) == 0.98*Kdc, w, 100))*1000;
s = tf('s');                                                                      %back to laplace s 
Q6.G = Kdc*(Q6.Wn)^2/(s^2+2*Q6.Z*Q6.Wn*s+(Q6.Wn)^2);                              %(V/V) xfer = Kdc*omega^2/(s^2+2*zeta*omega*s+omega^2)

%plots to check Q6
if 0
    [y t] = step(Q6.G);
    t = t*10^3;
    figure(1); clf;
    a2DSPlot(ABCDEFGH)  %want our approximation to be as close to the original plot as possible
    hold on
    plot(t,y);
    hold off
end

%heck all of the approximations until now on the same figure
if 0 
    [y t] = step(Q6.G);
    t = t*10^3;
    figure(1); clf;
    a2DSPlot(ABCDEFGH)
    hold on
    plot(t,y);
    [y t] = step(Q5.G);
    t = t*10^3;
    plot(t,y);
    [y t] = step(Q4.G);
    t = t*10^3;
    plot(t,y);
    [y t] = step(Q3.G);
    t = t*10^3;
    plot(t,y);
    [y t] = step(Q2.G);
    t = t*10^3;
    plot(t,y);
    hold off
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%Q7: Approximate #6

%plot to get the rise time
if 0
    [y t] = step(Q3.G);
    t = t*10^3;
    figure(1); clf;
    hold on
    plot(t,y);
    hold off
end

Tr2 = 0.21                                              %(ms) first time it gets to the final value before the peak
Q7.Tr = (Q1.Tr+0.21)/2;                                 %take their average
Q7.Te = Q7.Tr - Q1.Tr;                                  %error is the the difference between the average and the original value
Q7.Z = sqrt(log(os)^2/(pi^2 + log(os)^2));              %(pure) zeta = sqrt[(ln(overshoot)^2/(pi^2+ln(overshoot)^2)]
beta7 = sqrt(1-Q7.Z^2);                                 %(pure) beta = sqrt(1-zeta^2)
Q7.Wn = 1/(beta7*Q7.Tr*10^(-3))*(pi-atan(beta7/Q7.Z));  %(rad/s) omega with Tr = (pi-atan(beta/zeta))/(beta*Tr)
Q7.G = Kdc*(Q7.Wn^2)/(s^2+2*Q7.Z*Q7.Wn*s+Q7.Wn^2);      %(V/V) %(V/V) xfer = Kdc*omega^2/(s^2+2*zeta*omega*s+omega^2), Kdc was defined at the top

%plots to check Q7
if 0
    [y t] = step(Q7.G);
    t = t*10^3;
    figure(1); clf;
    a2DSPlot(ABCDEFGH)  %want our approximation to be as close to the original plot as possible
    hold on
    plot(t,y);
    hold off
end



