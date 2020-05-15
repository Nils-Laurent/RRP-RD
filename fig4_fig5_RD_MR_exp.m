close all;

%% signal definition
L = 4096;
t = (0:L-1)'/L;

B1 = log(2000);
phi1 = 200*t+2300*(t.^2)/2;
phi2 = 2000*t+exp(B1*t)/B1;
IF1 = 200+2300*t;
IF2 = 2000+exp(B1*t);
s1 = exp(2*1i*pi*phi1);
s2 = exp(2*1i*pi*phi2);
modes_exp_LC = transpose([s1, s2]);
IFs_exp_LC = transpose([IF1, IF2]);
sigma_exp_LC = 0.0241;

Nr = 2;
Nfft = 512;
cas = 1;
poly_degree = 5;

%% SNR (RD + MR) computation
clwin = 10;
SNR_IN = -10:5;
NRep = 10;

[SNR_NEW, SNR_C_RD, SNR_MB_RD, SNR_IF_NEW, SNR_IF_C_RD, SNR_IF_MB_RD] =...
    test_RD_MR(modes_exp_LC, IFs_exp_LC, clwin, sigma_exp_LC, Nfft, poly_degree, SNR_IN, NRep);

plot_SNR_modes(SNR_IN, SNR_NEW, SNR_C_RD, SNR_MB_RD);
plot_SNR_IFs(SNR_IN, SNR_IF_NEW, SNR_IF_C_RD, SNR_IF_MB_RD);