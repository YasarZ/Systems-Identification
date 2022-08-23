%% Collecting data from assigned flag 48 OE
flag = '48';
N = 5000;   %Number of samples
K = 10;     %Switching rate
A = 0.1;      %Ampltidue
u = A*dprbs(N,K);
[y, Ts] = process_OE(u, flag);

disp('*******************************');

%% Running Diagnostic on Flag
DC = -3;
y = y - DC;
u = vector (u, 'col');
y = vector (y, 'col');
dataid = iddata(y,u,Ts);

figure;
idplot(dataid(1:10));
title ('First 500 Samples of I/O data and check for saturation');

z = [y u];                          %input/output data
M = 500;                            %width of Hamming window
N = 1000;                           %length of data
Ts = 0.05                           %Sampling Time  
w = linspace(pi/(N/Ts),pi/Ts,N);    %Vector of frequency points
maxsize=[];
g = spa(z,M,w,maxsize,Ts);         %Frequency response information
bode(g);                           %Plotting the bode plot
legend ('SPA method');

% Test [2 3] and [3 3]

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

%Based on the Impulse Response, there are 8 delays in the system
%Based on the Impulse Response, we can except at least a second order system

%Hankel Test
figure;
hanktest(w, 8, 0);
figure;
hanktest(w, 8, 1);
legend("Noisy Data CRA");
%Hankel Test shows second order

%% Collecting data from assigned flag 48 OE
flag = '48';
N = 5000;   %Number of samples
K = 100;     %Switching rate
A = 0.8;      %Ampltidue
u = A*dprbs(N,K);
[y, Ts] = process_OE(u, flag);

disp('*******************************');
DC = -3;
y = y - DC;
u = vector (u, 'col');
y = vector (y, 'col');
dataid = iddata(y,u,Ts);
%% OE Model of Flag
datai = dataid(1:4700);
datav = dataid(4701:end);

nb = 2;
nf = 3;
nk = 8;
M_OE = oe(dataid, [nb nf nk]);

present (M_OE);

figure;
compare (dataid, M_OE);
disp('********************************************************');

% Fairly confident that the system is [2 3 8]



%% Validation for Flag
figure;
resid (M_OE, dataid);
residuals = resid(M_OE, dataid);
res = get(residuals, 'y');
res1 = residuals.y;
figure;
[S,x] = chisq(40,res);

yv = datav.y;
uv = datav.u;

figure;
crosscorr(res,uv,40);
[ccf] = crosscorr(res, uv, 40);

figure;
autocorr(res,40,1);
[acf2] = autocorr(res,40,1);

%% G(s) for Flag
Oc = d2c(M_OE);

a=[-1.8 -36;1 0];
b=[1;0];c=[ 25.74 51.48];d=0;
Tdel=0.35;Ts=0.05;
DC=-3;NoiseA=0.001;SAT=10;

[num, den] = ss2tf(a,b,c,d);
[numg, deng] = c2dm(num,den,Ts,'zoh');
Gd = tf(numg, deng, Ts, 'OutputDelay', 8);
Gs = d2c(Gd);

figure;
step(Gs,Gd,Oc,M_OE);
legend('CT Model','DT Model','OE Model in CT', 'OE Model');
figure;
bode(Gs,Gd,Oc,M_OE);
legend('CT Model','DT Model','OE Model in CT', 'OE Model');

num = [25.6 3202 6302];
den = [1 124.2 256.4 4405];
Gs = tf (num, den, 'IODelay', Ts);

test = zpk(Oc);
