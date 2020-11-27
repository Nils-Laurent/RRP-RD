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
% SNR_IN = [-8, -6, -4];
SNR_IN = -4;
NRep = 1;

% [SNR_LC] = test_RD_MR(modes_LC, IFs_LC, clwin, sigma_LC, Nfft, smooth_p, SNR_IN, NRep);
[SNR_cos] = test_RD_MR(modes_cos, IFs_cos, clwin, sigma_cos, Nfft, smooth_p, SNR_IN, NRep);
% [SNR_exp] = test_RD_MR(modes_exp, IFs_exp, clwin, sigma_exp, Nfft, smooth_p, SNR_IN, NRep);

R1_plot_fig6_7(SNR_IN, SNR_cos);