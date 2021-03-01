function [IF_FSST, m_FSST, FSST4, IF_new, m_simple, m_LCR, STFT_LCR, STFT] =...
    R1_GW_RRP(s_in, Fs, Nfft, sigma_L, sigma_Fs, smooth_p, N_Y)
L = length(s_in);
Nr = 1;

%% FSST4
% 1: verif meme fenetre

s1 = s_in;
% Te1=t1(2)-t1(1);
nv = log2(L);
if mod(nv,1)~=0
%     warning('The signal is not a power of two, zero padding to the next power');
    s1_pad = [s1.' zeros(1,2^(floor(log2(L))+1)-L)];
%     t1pad = max(t1)+Te1:Te1:max(t1)+(2^(floor(log2(L))+1)-L)*Te1;
%     t1_pad = [t1.' t1pad];
%     sw_pad = [sw zeros(1,2^(floor(log2(L))+1)-L)];
end
L_pad = length(s1_pad);
s_pad = s1_pad;
% s_pad = hilbert(s1_pad);

gamma = 10^(-3);
clwin = 10;
lambda = 0;
% L_pad = L;
% s_pad = s_noise;
ft =1:L_pad/32;bt=1:L_pad;

[STFT,FSST,FSST2,FSST3,FSST4,omega,omega2,omega3,tau2,tau3,phi22p,phi33p,phi44p] =...
    sstn(s_pad,gamma,sigma_L,ft,bt);
[Cs, Es] = exridge_mult_Noise(FSST4, Nr, lambda, clwin);
IF_FSST = (Cs(1, :) - 1)*Fs/Nfft;
d=10;
imf4 = sigma_L*recmodes(FSST4,Cs,d);
m_FSST = real(imf4);
% return;

%% RRP
% Fs1 = Fs;
% Fs = L;
sigma_s = sigma_Fs;

[g, Lh] = create_gaussian_window(Fs, Nfft, sigma_s);
% X_cmp = Lh+1:(L-Lh);
[STFT, omega, ~, QM, ~, tau] = FM_operators(s_in, Fs, Nfft, g, Lh, sigma_s);

% N_Y = 128;
STFT2 = STFT(1:N_Y, :);
QM2 = QM(1:N_Y, :);
omega2 = omega(1:N_Y, :);
tau2 = tau(1:N_Y, :);

% smooth_vec = [1 - 10^(-3), 1 - 10^(-6), 1 - 10^(-10)];
% smooth_p = 1 - 10^(-11);
% smooth_p = smooth_p;
% Ns = length(smooth_p);

% SNRs = zeros(Ns, 2);
% IFs_vec = zeros(Ns, L);

[Spl, ~] = R1_RRP_RD(STFT2, QM2, omega2, tau2, Fs, Nfft, Nr, sigma_s, smooth_p);

[m_simple, m_LCR, IF_new, STFT_LCR] = R1_MR_and_LCR_spl(STFT2, Spl, g, Lh, sigma_s, Nr, Nfft, Fs);
% IFs_vec(ind, :) = IF_new(1, :);
m_simple = real(m_simple);
m_LCR = real(m_LCR);
end

