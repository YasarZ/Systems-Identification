%% Collecting data from assigned flag 06 PEM
flag = '06';
N = 20000;   %Number of samples
K = 1;     %Switching rate
A = 2;      %Ampltidue
u = A*dprbs(N,K);
[y, Ts] = process_PEM(u, flag);

disp("*************************************************");

%% Running Diagnostic on Flag
DC = 4;
y = y - DC;
u = vector (u, 'col');
y = vector (y, 'col');
dataid = iddata(y,u,Ts);

figure;
idplot(dataid(1:500));
title ('First 500 Samples of I/O data and check for saturation');

z = [y u];                          %input/output data
M = 500;                            %width of Hamming window
N = 1000;                           %length of data
Ts = 0.05                           %Sampling Time
w = linspace(pi/(N/Ts),pi/Ts,N);    %Vector of frequency points
maxsize=[];
g2n = spa(z,M,w,maxsize,Ts);         %Frequency response information
bode(g2n);                           %Plotting the bode plot
legend ('SPA method');

%CRA 50 LAGS
L = 50;                             %number of lags
figure;
w = cra(dataid, L);
xlim([0 L/2]);
legend ('Data (CRA)');

%CRA 1000 LAGS
L = 1000;                           %number of lags
figure;
w = cra(dataid, L);
xlim([0 L/2]);
legend ('Data (CRA)');


figure;
hanktest(w, 7, 0);
figure;
hanktest(w, 7, 1);
legend("Noisy Data CRA");

%% OE MODEL
flag = '06';
N = 5000;   %Number of samples
K = 50;     %Switching rate
A = 1;      %Ampltidue
u = A*dprbs(N,K);
[y, Ts] = process_PEM(u, flag);
[y2, Ts2] = process_PEM(u,flag);

disp("*************************************************");
DC = 4;
y = y - DC;
u = vector (u, 'col');
y = vector (y, 'col');
y2 = y2 - DC;
y2 = vector (y2, 'col');
dataid = iddata(y,u,Ts);
dataid2 = iddata(y2,u,Ts);


datai = dataid(1:4700);
datav = dataid(4701:end);

nb = 2;
nf = 3;
nk = 7;
M_OE = oe(dataid, [nb nf nk]);

present (M_OE);

figure;
compare (dataid2, M_OE);
disp('********************************************************');
%% Validation
figure;
resid (M_OE, dataid);
residuals = resid(M_OE, dataid);
res = get(residuals, 'y');
figure;
[S,x] = chisq(40,res)

figure;
crosscorr(res,datav.u,100);
[ccf] = crosscorr(res, datav.u, 100);


figure;
[acf1] = acorrf(res, 100, 1);

figure;
[pacf1] = pacorrf(res, 100, 1);

figure;
parcorr(res,50);
%% Creating the Model
disp('----------------------Noise 1 Model--------------------');
Mbj = bj(datai, [2 0 1 3 7]);
%0 1 is good
present(Mbj);
figure;
compare(Mbj, datav, 1);


residuals = resid(Mbj, dataid);
res = get(residuals, 'y');


figure;
crosscorr(res,datav.u,40);
[ccf] = crosscorr(res, datav.u, 40);

figure;
[acf1] = acorrf(res, 100, 1);

figure;
[pacf1] = pacorrf(res, 100, 1);

%% Comparison 
Oc = d2c(Mbj);

a=[-4 -5 -2;1 0 0;0 1 0];
b=[1;0;0];c=[0 -10 10];d=0;
Tdel=0.3;
DC=4;NoiseA=0.1;SAT=15;
Ts=0.05;

[num, den] = ss2tf(a,b,c,d);
[numg, deng] = c2dm(num,den,Ts,'zoh');
Gd = tf(numg, deng, Ts, 'OutputDelay', 8);
Gs = d2c(Gd);

figure;
step(Gs,Gd,Oc,Mbj);
legend('CT Model','DT Model','OE Model in CT', 'OE Model');
figure;
bode(Gs,Gd,Oc,Mbj);
legend('CT Model','DT Model','OE Model in CT', 'OE Model');

num = [25.6 3202 6302];
den = [1 124.2 256.4 4405];
Gs = tf (num, den, 'IODelay', Ts);

%% G(s)
[A, B, C, D, F] = polydata(Mbj);
GD = tf(B,F,Ts);
GC = d2c(GD);
[nc, dc] = tfdata(GC, 'v');
Gmodel = tf(nc,dc,'IODelay', Ts);
test = zpk (GC);