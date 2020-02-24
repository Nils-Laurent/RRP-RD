close all

%% signal definition
L = 4096;
t = (0:L-1)'/L;

A = 45;
B = 4000;
phi = A*t+B*(t.^2)/2;
s_clean = exp(2*1i*pi*phi);

%%
sigma = 1/sqrt(B);

Nfft = 512;
cas = 1;

noise = randn(L,1)+1i*randn(L,1);
s_noise = sigmerge(s_clean, noise, -10);

Wnoise = s_noise - s_clean;
std(real(Wnoise))
std(imag(Wnoise))
pause

[g, Lg] = create_gaussian_window(L, Nfft, sigma);

%% 2nd order computation
[TFR_noise, omega, omega2, q] = FM_operators(s_noise, Nfft, g, Lg, sigma);

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
