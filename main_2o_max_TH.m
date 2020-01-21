close all

%% signal definition
L = 4096;
t = (0:L-1)'/L;

A = 45;
B = 4000;
phi = A*t+B*(t.^2)/2;
s_clean = (3/4)*exp(2*1i*pi*phi);
phip = A + B*t;
phipp = B*ones(L, 1);

%%
sigma = 1/sqrt(B);

Nfft = 512;
cas = 1;

noise = randn(L,1)+1i*randn(L,1);
s_noise = sigmerge(s_clean, noise, -10);
[g, Lg] = create_gaussian_window(L, Nfft, sigma);

%% 2nd order computation
[TFR_noise, omega, omega2, q] = FM_operators(s_noise, Nfft, g, Lg, sigma);
exridge_new(TFR_noise, sigma, q, 2);