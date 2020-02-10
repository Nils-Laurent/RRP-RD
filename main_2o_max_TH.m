close all

%% signal definition
L = 4096;
t = (0:L-1)'/L;

A = 45;
B = 4000;
phi = A*t+B*(t.^2)/2;
s_clean = (3/4)*exp(2*1i*pi*phi);

%%
sigma = 1/sqrt(B);

Nfft = 512;
cas = 1;

N_noise = 100;
R_set = 8:16:1024;
NR = length(R_set);
E = zeros(N_noise, NR);
for n=1:N_noise
    fprintf("n = %d/%d\n", n, N_noise);
    noise = randn(L,1)+1i*randn(L,1);
    s_noise = sigmerge(s_clean, noise, -10);

    [g, Lg] = create_gaussian_window(L, Nfft, sigma);

    %% 2nd order computation
    [TFR_noise, omega, omega2, q] = FM_operators(s_noise, Nfft, g, Lg, sigma);

    [C_test, E_curve] = exridge_n2(TFR_noise, q, 2, R_set);
    E(n, :) = E_curve;
end

E_smooth = mean(E, 1);
plot(R_set, E_smooth);
