%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% MATLAB QUIZ - ELEC 341 - IDIL BIL
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
%%%Q1: RLC Circuit

% plot xmtPlot(ABCDEFGH) on the command line and find the second order approximation times 
Tr = 48.2253 * 10^(-3);         %(s) first time it gets to the final value before the peak
Tp = 76.5932 * 10^(-3);         %(s) first time it gets to the highest value
%no need for settle because you won't use it in your calculations
Pos = (1.28298-1)/1 * 100;      %(percent) (Tpeak y - Tsettle y)/Tsettle y * 100
FV = 1.0;                       %gain = y value of Tsettle
overshoot = Pos/100;

Zeta = sqrt((log((overshoot)))^2/(pi^2+(log((overshoot)))^2));      %(pure) zeta = sqrt[(ln(overshoot)^2/(pi^2+ln(overshoot)^2)]
Beta = sqrt(1-Zeta^2);                                              %(pure) beta = sqrt(1-zeta^2)
wn = (1/(Beta*Tr)) * (pi-atan(Beta/Zeta));                          %(rad/s) omega with Trise = (pi-atan(beta/zeta))/(beta*Tr)
s = tf('s');                                                        %define s for laplace
tf = FV*(wn^2)/(s^2+2*Zeta*wn*s+wn^2);                              %(V/V) xfer = Kdc*omega^2/(s^2+2*zeta*omega*s+omega^2)

if 0
    figure(1);
    clf;
    grid on;
    xmtPlot(ABCDEFGH);              %check if the plots overlap if not try with another time like Tpeak
    [test_y test_t] = step(tf);     %[y x] = step(xfer)
    test_t = test_t*10^3;           %change timescale to ms because thats what xmtPlot uses
    hold on;
    plot(test_t, test_y, 'r--' ,'LineWidth', 3);
end

%write a voltage dividor function for C (thats how V is defined)
%V = [1/sC]/[R+(1/sC)+sL] -> expand with s -> [1/C]/[s^2L+sR+1/C] -> expand with (1/L) to match tf -> [1/(LC)]/[s^2+R/Ls+1/(LC)]
%look at tf from the command line
%compute C and L
R = 50;             %(ohms)
L = R/32.54;        %(henry)
C = 1/(L*1905);     %(farad)

Q1.L = L;           %(henry)
Q1.C = C*10^(6);    %(microfarad)

%write the tf with the values found
test = (1/(Q1.L*Q1.C*10^(-6)))/(s^2+R/Q1.L*s+(1/(Q1.L*Q1.C*10^(-6))));     %tf = [1/(LC)]/[s^2+R/Ls+1/(LC)]

if 0
    [y t] = step(test);
    t = t*10^3;                     %change timescale to ms because thats what xmtPlot uses   
    plot(t,y,'b-', 'LineWidth', 2); %check if the plots overlap if not your impedances are wrong
    hold off;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%Q2: State-Space Matrixes

%variable declaration for generator spesifications
Rw = A/5;       %(ohms)
Lw = (B+C)/10;  %(henry)
Jg = (D+E)/6;   %(Nms^2/rad)
Bg = (F+G)/5;   %(Nms/rad)
Km = (H)/2;     %(Nm/A)

%[sv; si; sw] = A * [v; i; w] + B * [Ta; ](inputs)
matrix_A = [0 1/(Q1.C*10^(-6)) 0;
            -1/(L+Lw) -(Rw+R)/(L+Lw) Km/(L+Lw);
            0 -Km/Jg -Bg/Jg];
matrix_B = [0 0 1/Jg].';

Q2.A = matrix_A;
Q2.B = matrix_B;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%Q3: Transfer Function

%[V; i](outputs) = C * [V; i; w] + D * [Ta; ](inputs)
%for D -> output number = number of rows, input number = number of columns
matrix_C = [1 0 0;
            0 1 0];
matrix_D = [0; 0];

ID_matrix = eye(3);
temp = s*ID_matrix-matrix_A;
STM = inv(temp);                              %system transition matrix = (s*identity matrix - A)^(-1)
TF_matrix = matrix_C*STM*matrix_B+matrix_D;   %transfer function = C*STM*B+D

Q3.Xv = TF_matrix(1);
Q3.Xi = TF_matrix(2);

if 
    figure(3); clf;
    step(Q3.Xv);        %check if its stable

    figure(4); clf;
    step(Q3.Xi);        %check if its stable
end
