close all;

addpath('./RRP_alg/');
addpath('./test/');

%% global var

L = 4096;
t = (0:L-1)'/L;

%% Signal def
% B = 3*L/4;
% phi_LC = L/8*t+B*(t.^2)/2;
% s_LC = exp(2*1i*pi*phi_LC);
% sigma_LC = 1/sqrt(B);
% 
% B = floor(2*L/(2*pi));
% phi_cos = L/2*t+B/(2*pi)*cos(2*pi*t);
% s_cos = exp(2*1i*pi*phi_cos);
% sigma_cos = 0.0142;
% 
% B = log(4096 - 510);
% phi_exp = 500*t+exp(B*t)/B;
% s_exp = exp(2*1i*pi*phi_exp);
% sigma_exp = 0.028;

phi1 = 256*t+2700*(t.^2)/2;
phi2 = 768*t+3000*(t.^2)/2;
s1 = exp(2*1i*pi*phi1);
s2 = exp(2*1i*pi*phi2);
s_LC = s1 + s2;
sigma_LC = 0.0188;

phi1 = 1400*t+1350/(2*pi)*cos(2*pi*t + pi/2);
phi2 = 3400*t+550/(2*pi)*cos(2*pi*t);
s1 = exp(2*1i*pi*phi1);
s2 = exp(2*1i*pi*phi2);
s_cos = s1 + s2;
sigma_cos = 0.0175;

B1 = log(2000);
phi1 = 200*t+2300*(t.^2)/2;
phi2 = 2000*t+exp(B1*t)/B1;
s1 = exp(2*1i*pi*phi1);
s2 = exp(2*1i*pi*phi2);
s_exp = s1 + s2;
sigma_exp = 0.0241;

%% test
SNR_in = -10;
Nfft = 512;
smooth_p = 1 - 10^(-4);
Nr = 2;

% 3 FIG sig. MCS
% 3 Courbes (spline, std sup, std inf)

%% LC
% noise = randn(L, 1)+1i*randn(L, 1);
load('mat/noise_R1_TFR_RD_LC.mat');
s_noise = sigmerge(s_LC, noise, SNR_in);
[g, Lh] = create_gaussian_window(L, Nfft, sigma_LC);
[STFT, omega, ~, QM, ~, tau] = FM_operators(s_noise, L, Nfft, g, Lh, sigma_LC);
[Spl_LC, ~] = R1_RRP_RD(STFT, QM, omega, tau, L, Nfft, Nr, sigma_LC, smooth_p);

Spl = Spl_LC;
sigma_s = sigma_LC;
IF1 = fnval(Spl(1).spline, t);
DF1 = fnval(fnder(Spl(1).spline), t);
R1 = 1/(sqrt(2*pi)*sigma_s)*sqrt(1 + sigma_s^4*DF1.^2);
IF2 = fnval(Spl(2).spline, t);
DF2 = fnval(fnder(Spl(2).spline), t);
R2 = 1/(sqrt(2*pi)*sigma_s)*sqrt(1 + sigma_s^4*DF2.^2);

% save('noise_R1_TFR_RD_LC.mat', 'noise');
fname = 'fig_R1_TFR_RD_LC';
R1_plot_fig5(STFT, IF1, R1, IF2, R2, fname);
% return;

%% cos
% noise = randn(L, 1)+1i*randn(L, 1);
load('mat/noise_R1_TFR_RD_cos.mat');
s_noise = sigmerge(s_cos, noise, SNR_in);
[g, Lh] = create_gaussian_window(L, Nfft, sigma_cos);
[STFT, omega, ~, QM, ~, tau] = FM_operators(s_noise, L, Nfft, g, Lh, sigma_cos);
[Spl_cos, ~] = R1_RRP_RD(STFT, QM, omega, tau, L, Nfft, Nr, sigma_cos, smooth_p);

Spl = Spl_cos;
sigma_s = sigma_cos;
IF1 = fnval(Spl(1).spline, t);
DF1 = fnval(fnder(Spl(1).spline), t);
R1 = 1/(sqrt(2*pi)*sigma_s)*sqrt(1 + sigma_s^4*DF1.^2);
IF2 = fnval(Spl(2).spline, t);
DF2 = fnval(fnder(Spl(2).spline), t);
R2 = 1/(sqrt(2*pi)*sigma_s)*sqrt(1 + sigma_s^4*DF2.^2);

% save('noise_R1_TFR_RD_cos.mat', 'noise');
fname = 'fig_R1_TFR_RD_cos';
R1_plot_fig5(STFT, IF1, R1, IF2, R2, fname);
% return;

%% exp
% noise = randn(L, 1)+1i*randn(L, 1);
load('mat/noise_R1_TFR_RD_exp.mat');
s_noise = sigmerge(s_exp, noise, SNR_in);
[g, Lh] = create_gaussian_window(L, Nfft, sigma_exp);
[STFT, omega, ~, QM, ~, tau] = FM_operators(s_noise, L, Nfft, g, Lh, sigma_exp);
[Spl_exp, ~] = R1_RRP_RD(STFT, QM, omega, tau, L, Nfft, Nr, sigma_exp, smooth_p);

Spl = Spl_exp;
sigma_s = sigma_exp;
IF1 = fnval(Spl(1).spline, t);
DF1 = fnval(fnder(Spl(1).spline), t);
R1 = 1/(sqrt(2*pi)*sigma_s)*sqrt(1 + sigma_s^4*DF1.^2);
IF2 = fnval(Spl(2).spline, t);
DF2 = fnval(fnder(Spl(2).spline), t);
R2 = 1/(sqrt(2*pi)*sigma_s)*sqrt(1 + sigma_s^4*DF2.^2);

% save('noise_R1_TFR_RD_exp.mat', 'noise');
fname = 'fig_R1_TFR_RD_exp';
R1_plot_fig5(STFT, IF1, R1, IF2, R2, fname);