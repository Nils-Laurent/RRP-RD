close all;

%% signal definition
L = 4096;
t = (0:L-1)'/L;

% B = 3*L/4;
% phi_s = L/8*t+B*(t.^2)/2;
% s_in = exp(2*1i*pi*phi_s);
% sigma_s = 1/sqrt(B);

B = floor(2*L/(2*pi));
phi_s = L/2*t+B/(2*pi)*cos(2*pi*t);
s_in = exp(2*1i*pi*phi_s);
sigma_s = 0.0142;
% 
% phip_s = L/2-B*sin(2*pi*t);
% phipp_s = -2*pi*cos(2*pi*t);
% 
% std_s = 1/(sqrt(2*pi)*sigma_s)*sqrt(1 + sigma_s^4*phipp_s.^2);

%% STFT
Nfft = 512;
Nr = 1;
cas = 1;

fx = (0:Nfft-1)*L/Nfft;

noise = randn(L,1)+1i*randn(L,1);
sn_LC = sigmerge(s_in, noise, -10);
[g, Lh] = create_gaussian_window(L, Nfft, sigma_s);

[STFT, omega, ~, QM, Vxgp, tau] = FM_operators(sn_LC, Nfft, g, Lh, sigma_s);

%% RRP
degree_WPF = 5;
R1_RRP_RD(STFT, QM, tau, omega, sigma_s, Nr, degree_WPF);