close all

%% signal definition
L = 4096;
t = (0:L-1)'/L;

f2_amp = 3000;
phi1 = 200*t+f2_amp*(t.^2)/2;
phi2 = 265*t+f2_amp*(t.^2)/2;

s_clean = exp(2*1i*pi*phi1) + exp(2*1i*pi*phi2);
NRidges = 2;

Nr = 2;

%%
sigma = 1/sqrt(f2_amp);
% eta_lim = 1/sqrt(2*pi)*sqrt(1/sigma^2 + sigma^2*phipp^2);

Nfft = 512;
cas = 1;

noise = (2.237)*(randn(L,1)+1i*randn(L,1));
%noise = (0.001)*(randn(L,1)+1i*randn(L,1));
% s_noise = sigmerge(s_clean, noise, -5);
s_noise = s_clean + noise;

[g, Lg] = create_gaussian_window(L, Nfft, sigma);

%% 2nd order computation
[TFR_noise, omega, omega2, q] = FM_operators(s_noise, Nfft, g, Lg, sigma);

% figure;
% imagesc(1:L, 1:Nfft, abs(TFR_noise));
% set(gca,'ydir','normal');
% axis square
% pause

C = 2;
[Cs] = exridge_n2_MCS(TFR_noise, q, C, Nr, sigma);

Cs(Cs == 0) = nan;

figure;
imagesc(1:L, 1:Nfft, abs(TFR_noise));
set(gca,'ydir','normal');
axis square
hold on;
plot(1:L, Cs);
hold off;
