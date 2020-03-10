close all

%% signal definition
L = 4096;
t = (0:L-1)'/L;

A = 45;
B = 4000;
phi = A*t+B*(t.^2)/2;
s_clean = exp(2*1i*pi*phi);

%%
sigma_s = 1/sqrt(B);

Nfft = 512;
cas = 1;

noise = randn(L,1)+1i*randn(L,1);
s_noise = sigmerge(s_clean, noise, -10);

% Wnoise = s_noise - s_clean;
% std(real(Wnoise))
% std(imag(Wnoise))
% pause

[g, Lg] = create_gaussian_window(L, Nfft, sigma_s);

%% 2nd order computation
[TFR_noise, omega, omega2, q] = FM_operators(s_noise, Nfft, g, Lg, sigma_s);

[Cr] = exridge_n2(TFR_noise, q, sigma_s);

Cr(Cr == 0) = nan;

figure;
imagesc(1:L, 1:Nfft, abs(TFR_noise));
set(gca,'ydir','normal');
colormap(flipud(gray));
axis square
hold on;
plot(Cr, 'r');
hold off;
