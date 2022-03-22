addpath('./toolbox/TF_Toolbox/');
addpath('./RRP_alg/');

close all;

%% signal samples
Lx = 1024;
Tx = (0:(Lx-1))/Lx;

%% signal definition
B = 400;
phi1 = 100*Tx + B*Tx.^2;
s1 = exp(2i*pi*phi1);

%% add complex noise
noise = randn(1, Lx) + 1i*randn(1, Lx); % gaussian white noise

snr_in = -10;
sn = add_noise(s1, noise, snr_in);

%% time frequency parameters
sigma_opt = 1/sqrt(B); % minimizing renyi entropy
[g, ~] = gauss_win(Lx, sigma_opt);

Nfft = 1024;
Fx = (0:(Nfft-1))*Lx/Nfft;

%% second order synchrosqueezing
[STFT_sn, TFR] = sst2(sn, sigma_opt, Nfft);

figure;
imagesc(Tx, Fx, abs(STFT_sn));
set(gca,'YDir','normal');
xlabel("time");
ylabel("frequency");
title("short time Fourier transform");

%% ridge detection

% P is the smoothness parameter for the weigthed spline approximation
% typicall value for synthetic signals at -10dB is 10^(-4)
P = 10^(-4);
[Spl, ~] = RRP_RD(STFT_sn, TFR.q_hat, TFR.omega1_hat, TFR.tau, P);

figure;
imagesc(Tx, Fx, abs(STFT_sn));
set(gca,'YDir','normal');
xlabel("time");
ylabel("frequency");
hold on;
plot(Tx, ppval(Spl(1).spline, Tx), 'r-');
hold off;
title("result of the ridge detection");