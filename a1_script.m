%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% MATLAB ASSIGNMENT 1 - ELEC 341 - IDIL BIL
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
%%%Q1: CCT Analysis

%variables for the voltage amplifier circuit
R1 = 10*A;      %(Ohms)
L1 = B*10^(-3); %(H)
C1 = C*10^(-6); %(F)
R2 = 10*D;      %(Ohms)
L2 = E*10^(-3); %(H)
C2 = F*10^(-6); %(F)
Vin = 18;       %(V)

%impedances
%R -> R, L -> sL, C = 1/(sC)
s = tf('s');     %define s for laplace
Z_R1 = R1;
Z_R2 = R2;
Z_C1 = 1/(s*C1);
Z_C2 = 1/(s*C2);
Z_L1 = s*L1;
Z_L2 = s*L2;

%transfer functions
%op amp 1
xfer1_num = 1/[(1/Z_R1) + (1/Z_L1) + (1/Z_C1)]; %resistor, inductor and capacitor in parallel
xfer1_denom = 1/[(1/Z_R1) + (1/Z_C1)];          %resistor and capacitor in parallel
xfer1 = minreal(-xfer1_num/xfer1_denom);        %xfer = -Z2/Z1   
%op amp 2
xfer2_num = 1/[(1/Z_R2) + (1/Z_C2)];            %resistor and capacitor in parallel
xfer2_denom = 1/[(1/Z_R2) + (1/Z_L2)];          %resistor and inductor in parallel
xfer2 = minreal(-xfer2_num/xfer2_denom);        %xfer = -Z2/Z1  

xfer_t = minreal(xfer1*xfer2);  %cascaded xfers = xfer1*xfer2
Q1.G = xfer_t;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%Q2: Impulse & Envelope

scaled_tf = Q1.G * Vin;     %scale with Vin because Vin = 1 for transfer functions
figure(1);
clf;
step(scaled_tf);            %compute and plot step response
grid on;

%longer way 
% [ys ts] = step(xfer_t);             %[y x] = step(x), compute the step function of the overall transfer function
% ys = ys*Vin;                        %scale with Vin because Vin = 1 for transfer functions
% figure(1);
% clf;
% xfer_t_step = plot(ts, ys, 'b-');   %plot(x, y, blue line graph) the step response
% set(xfer_t_step, 'LineWidth', 2);   %set line width
% grid on;
% title('Step Response');             %edit the title
% xlabel('Time (sec)');               %edit the xlabel
% ylabel('OP Voltage (V)');           %edit the ylabel

%computer the inverse laplace transform to find step response in time domain
step_scaled_tf_s = scaled_tf * (1/s);   %compute the step function
[num den] = tfdata(step_scaled_tf_s);
syms s t;                               %turns s and t to letters/symbols (to define a value or do ilaplace)
snum = poly2sym(num,s);
sden = poly2sym(den,s);
step_scaled_sym_s = snum/sden;
step_scaled_sym_t = ilaplace(step_scaled_sym_s);
step_scaled_sym_t = vpa(step_scaled_sym_t);

%this part is adjusted with the values calculated in the script
if 0
    t = 0:0.000050:0.1
    %ilaplace = X * e^(-A*t) * cos(w*t) + Y * e^(-A*t) * sin(w*t)
    %test = step_scaled_sym_t
    test = 0.20021628302178175878274398403899*exp(-649.35064935064932792961997279708*t) - 16.700216283021795619981169277017*exp(-320.51282051282042177515137325136*t).*cos(2624.9340039120802530566631525966*t) + 2.4631814849519443149128549559276*exp(-320.51282051282042177515137325136*t).*sin(2624.9340039120802530566631525966*t) + 16.500000000000013861198425292978;
    figure(7); 
    clf;
    plot(test);                         %needs to be enveloped
    % env = FV(the number at the end) +/- K*e^(-A*t)
    test_env1 = 16.500000000000013861198425292978 + sqrt(16.700216283021795619981169277017.^2+2.4631814849519443149128549559276^2)*exp(-320.51282051282042177515137325136*t);
    test_env2 = 16.500000000000013861198425292978 - sqrt(16.700216283021795619981169277017.^2+2.4631814849519443149128549559276^2)*exp(-320.51282051282042177515137325136*t);
    hold on;
    plot(test_env1, 'Color', 'r');      %make sure the envelopes fit
    plot(test_env2, 'Color', 'r');      %make sure the envelopes fit
end

%K = sqrt(X^2 + Y^2)
Q2.K = sqrt(16.700216283021795619981169277017.^2+2.4631814849519443149128549559276^2);
%coefficient of t
Q2.A = 320.51282051282042177515137325136;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%Q3: Approximate Impulse
%couldn't do it
if 0
    figure(2);
    clf;
    impulse(Q1.G);
    grid on;
end

Q3.Tl = 0;
Q3.Vh = 0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%Q4: CCT Analysis
Iin = H;

%after nodal analysis and hp prime (takes too long)
%use the matrix method thought in class instead
% I(input sources) = Y(admittances) * V(nodes) -> V = inv(Y) * I -> V(1) = tf of V1
% inputs -> 1 fot transfer functions (reccomended for single input), inputs -> X for scaled (reccomended for multiple inputs), inputs -> X/s for step responses
q4_num = [C1*C2^2*L1*L2^2 C1*C2^2*L2^2*R1+C1*C2^2*L1*L2*R2 C1*C2^2*L2*R1*R2+C2^2*L2^2+3*C1*C2*L1*L2 C2^2*L2*R2+2*C1*C2*L1*R2+3*C1*C2*L2*R1 2*C1*C2*R1*R2+C1*L1+3*C2*L2 C1*R1+2*C2*R2 1];
q4_den = [C1*C2^2*L1*L2^2 C1*C2^2*L2^2*R1+2*C1*C2^2*L1*L2*R1+2*C1*C2^2*L1*L2*R2 2*C1*C2^2*L1*R1*R2+2*C1*C2^2*L2*R1*R2+C2^2*L2^2+2*C2^2*L1*L2+3*C1*C2*L1*L2 2*C2^2*L1*R2+2*C2^2*L2*R2+C1*C2*L1*R1+2*C1*C2*L1*R2+3*C1*C2*L2*R1 2*C1*C2*R1*R2+C1*L1+C2*L1+3*C2*L2 C1*R1+2*C2*R2 1];
q4_tf = tf(q4_num,q4_den);

Q4.G = minreal(q4_tf);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%Q5: CCT Step & Envelope
%same as Q2

%V(4)/Z10 = Iout -> transfer functions

%Step response plot
[ys ts] = step(q4_tf);
figure(5);
clf;
q4_tf_t_step = plot(ts, ys, 'b-');
set(q4_tf_t_step, 'LineWidth', 2);
grid on;
title('Step Response');
xlabel('Time (sec)');
ylabel('OP Voltage (V)');

%Computer the inverse laplace transform to find step response in time domain
s = tf('s');
step_tf_4_s = q4_tf * (1/s);
[num den] = tfdata(step_tf_4_s);
syms s t;
snum = poly2sym(num,s);
sden = poly2sym(den,s);
step_tf_4_sym_s = snum/sden;
step_tf_4_sym_t = ilaplace(step_tf_4_sym_s);
step_tf_4_sym_t = vpa(step_tf_4_sym_t);

%this part is adjusted with the values calculated in the script
if 0
    t=0:0.000050:0.1;
    %ilaplace = X * e^(-A*t) * cos(w*t) + Y * e^(-A*t) * sin(w*t)
    %test4 = step_tf_4_sym_t
    test4 = 0.7406307137833127547829480944764020570004686656254386146089210544*exp(-36931.25714265215662944856682936766383200329376326058214176590869*t) + 0.0007556907215699452188897579225468720005929194303456591140432114182*exp(-9867.2975573616046561759301493223392158776925624566125271532128*t) - 0.009020322003780758564390661898267962592953047789548429290170293567*exp(-335.8604800108101200472541678726918817675154915799357082060100675*t) - 0.003078588150805864585788300293504489878758178498287010575428136781*exp(-681.4864207860603469341463789818323027378607734392934034009343027*t) - 0.7292874943502960768516588902071764765293503587679488338573658355*exp(-118.0232255687083731307337836548790595368081919933535483021610205*t).*cos(1876.505095682747303894905234159288463189668839318812150997542543*t) + 0.06705144205192394579485781169117633885631722760484908173050962833*exp(-118.0232255687083731307337836548790595368081919933535483021610205*t).*sin(1876.505095682747303894905234159288463189668839318812150997542543*t) + 1.0;
    figure(6); 
    clf;
    plot(test4);                        %needs to be enveloped
    % env = FV(the number at the end) +/- K*e^(-A*t)
    test_env1_4 = 1 + sqrt(0.06705144205192394579485781169117633885631722760484908173050962833^2+0.7292874943502960768516588902071764765293503587679488338573658355^2)*exp(-118.0232255687083731307337836548790595368081919933535483021610205*t);
    test_env2_4 = 1 - sqrt(0.06705144205192394579485781169117633885631722760484908173050962833^2+0.7292874943502960768516588902071764765293503587679488338573658355^2)*exp(-118.0232255687083731307337836548790595368081919933535483021610205*t);
    hold on;
    plot(test_env1_4, 'Color', 'r');    %make sure the envelopes fit
    plot(test_env2_4, 'Color', 'r');    %make sure the envelopes fit
end

%K = sqrt(X^2 + Y^2)
Q5.K = sqrt(0.06705144205192394579485781169117633885631722760484908173050962833^2+0.7292874943502960768516588902071764765293503587679488338573658355^2);
%coefficient of t
Q5.A = 118.0232255687083731307337836548790595368081919933535483021610205;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%Q6: Settle Time

%use stepinfo(function) to get all the times (s)
%stepinfo(Q4.G)
Q6.Ts = 0.0335;
