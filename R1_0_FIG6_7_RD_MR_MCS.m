close all;

%% global var
L = 4096;
t = (0:L-1)'/L;

%% Linear Chirp MCS

phi1 = 256*t+2700*(t.^2)/2;
phi2 = 768*t+3000*(t.^2)/2;
IF1 = 256+2700*t;
IF2 = 768+3000*t;
s1 = exp(2*1i*pi*phi1);
s2 = exp(2*1i*pi*phi2);
modes_LC = [s1, s2].';
IFs_LC = [IF1, IF2].';
sigma_LC = 0.0188;

%% Cosine MCS

phi1 = 1400*t+1350/(2*pi)*cos(2*pi*t + pi/2);
phi2 = 3400*t+550/(2*pi)*cos(2*pi*t);
if1 = 1400-1350*sin(2*pi*t + pi/2);
if2 = 3400-550*sin(2*pi*t);
s1 = exp(2*1i*pi*phi1);
s2 = exp(2*1i*pi*phi2);
modes_cos = [s1, s2].';
IFs_cos = [if1, if2].';
sigma_cos = 0.0175;

%% Exp MCS

B1 = log(2000);
phi1 = 200*t+2300*(t.^2)/2;
phi2 = 2000*t+exp(B1*t)/B1;
IF1 = 200+2300*t;
IF2 = 2000+exp(B1*t);
s1 = exp(2*1i*pi*phi1);
s2 = exp(2*1i*pi*phi2);
modes_exp = [s1, s2].';
IFs_exp = [IF1, IF2].';
sigma_exp = 0.0241;

%% Test
clwin = 10;
Nfft = 512;
smooth_p = [1 - 10^(-2), 1 - 10^(-3), 1 - 10^(-4)];
SNR_IN = [-8, -6, -4];
NRep = 1;

% [SNR_LC] = test_RD_MR(modes_LC, IFs_LC, clwin, sigma_LC, Nfft, smooth_p, SNR_IN, NRep);
% [SNR_cos] = test_RD_MR(modes_cos, IFs_cos, clwin, sigma_cos, Nfft, smooth_p, SNR_IN, NRep);
% [SNR_exp] = test_RD_MR(modes_exp, IFs_exp, clwin, sigma_exp, Nfft, smooth_p, SNR_IN, NRep);

%% fig var
B1 = [0 0.4470 0.7410];
B3 = [0.3010 0.7450 0.9330];
B2 = (B3 + B1)/2;

R_MB = [0.6350 0.0780 0.1840];
Y_Cl = [0.9290 0.6940 0.1250];

%% fig 6 : RD
% Deux lignes
% 6 FIG : (f1, f2), (Signal x3 MCS)
% 5 courbes : (p = 10^(-2, -3, -4), MB RD, S RD)

figure;
hold on;
plot(SNR_IN, SNR_cos.RD.Cl(1, :), ':',...
    'Color', Y_Cl);
plot(SNR_IN, SNR_cos.RD.MB(1, :), '-',...
    'Color', R_MB);
plot(SNR_IN, SNR_cos.RD.New(1, :, 1), '-*',...
    'Color', B1);
plot(SNR_IN, SNR_cos.RD.New(1, :, 2), '-s',...
    'Color', B2);
plot(SNR_IN, SNR_cos.RD.New(1, :, 3), '-o',...
    'Color', B3);
hold off;

%% fig 7 : MR LCR
% Deux lignes

% 3 FIG : MR 3 Sig MCS
% 8 Courbes : (f1, f2), (p_min, p_max, MB-MR, SR-MR)

figure;
hold on;
plot(SNR_IN, SNR_cos.MR.Cl(1, :), ':',...
    'Color', Y_Cl);
plot(SNR_IN, SNR_cos.MR.MB(1, :), '-',...
    'Color', R_MB);
plot(SNR_IN, SNR_cos.MR.New(1, :, 1), '-*',...
    'Color', B1);
plot(SNR_IN, SNR_cos.MR.New(1, :, 3), '-o',...
    'Color', B3);
plot(SNR_IN, SNR_cos.MR.Cl(2, :), '-.',...
    'Color', Y_Cl);
plot(SNR_IN, SNR_cos.MR.MB(2, :), '--',...
    'Color', R_MB);
plot(SNR_IN, SNR_cos.MR.New(2, :, 1), '--*',...
    'Color', B1);
plot(SNR_IN, SNR_cos.MR.New(2, :, 3), '--o',...
    'Color', B3);
hold off;

% 3 FIG : LCR 3 Sig MCS
% 7 Courbes : (f1, f2), (smooth_p x3) + SR-MR-LCR

figure;
hold on;
plot(SNR_IN, SNR_cos.LCR.Cl(1, :), ':',...
    'Color', Y_Cl);
plot(SNR_IN, SNR_cos.LCR.New(1, :, 1), '-*',...
    'Color', B1);
plot(SNR_IN, SNR_cos.LCR.New(1, :, 2), '-s',...
    'Color', B2);
plot(SNR_IN, SNR_cos.LCR.New(1, :, 3), '-o',...
    'Color', B3);
plot(SNR_IN, SNR_cos.LCR.Cl(2, :), '-.',...
    'Color', Y_Cl);
plot(SNR_IN, SNR_cos.LCR.New(2, :, 1), '--*',...
    'Color', B1);
plot(SNR_IN, SNR_cos.LCR.New(2, :, 2), '--s',...
    'Color', B2);
plot(SNR_IN, SNR_cos.LCR.New(2, :, 3), '--o',...
    'Color', B3);
hold off;