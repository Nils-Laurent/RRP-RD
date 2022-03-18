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
    sstn_test(s_pad,gamma,sigma_L,ft,bt);
[Cs, Es] = exridge_mult_Noise(FSST4, Nr, lambda, clwin);
IF_FSST = (Cs(1, :) - 1)*Fs/Nfft;
d=10;
imf4 = sigma_L*recmodes(FSST4,Cs,d);
m_FSST = real(imf4);
% return;

%% RRP

[STFT, TFR] = sst2(s_in, sigma_L, Nfft);
QM = TFR.q_hat;
omega = TFR.omega1_hat;
tau = TFR.tau;

% N_Y = 128;
STFT2 = STFT(1:N_Y, :);
QM2 = QM(1:N_Y, :);
omega2 = omega(1:N_Y, :);
tau2 = tau(1:N_Y, :);

[Spl, ~] = RRP_RD(STFT2, QM2, omega2, tau2, smooth_p, Nr, 'Nfft', Nfft);

[g, Lh] = gauss_win(L, sigma_L);
[m_simple, m_LCR, IF_new, STFT_LCR] =...
    R1_MR_and_LCR_spl(STFT2, Spl, g, Lh, sigma_L, Nr, Nfft, L);

m_simple = real(m_simple);
m_LCR = real(m_LCR);
IF_new = Fs/L*IF_new;
end

