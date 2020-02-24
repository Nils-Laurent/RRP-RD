close all

%% signal definition
L = 4096;
t = (0:L-1)'/L;

A = 45;
B = 4000;
phi1 = A*t+3/4*B*(t.^2)/2;
phi2 = 10*A*t+4/5*B*(t.^2)/2;
s_clean = exp(2*1i*pi*phi1) + exp(2*1i*pi*phi2);

%%
sigma = 1/sqrt(B);
% eta_lim = 1/sqrt(2*pi)*sqrt(1/sigma^2 + sigma^2*phipp^2);

Nfft = 512;
cas = 1;

noise = (2.237)*(randn(L,1)+1i*randn(L,1));
%s_noise = sigmerge(s_clean, noise, -5);
s_noise = s_clean + noise;

[g, Lg] = create_gaussian_window(L, Nfft, sigma);

%% 2nd order computation
[TFR_noise, omega, omega2, q] = FM_operators(s_noise, Nfft, g, Lg, sigma);

figure;
imagesc(1:L, 1:Nfft, abs(TFR_noise));
set(gca,'ydir','normal');
axis square
pause

[C_test, curves, energies] = exridge_n2(TFR_noise, q, 2);

C_test(C_test == 0) = nan;
curves(curves == 0) = nan;

figure;
plot(energies);

figure;
imagesc(1:L, 1:Nfft, abs(TFR_noise));
set(gca,'ydir','normal');
axis square
hold on;
plot(C_test, 'r');
hold off;

figure;
imagesc(1:L, 1:Nfft, abs(TFR_noise));
set(gca,'ydir','normal');
axis square
hold on;
plot(1:L, curves);
hold off;
