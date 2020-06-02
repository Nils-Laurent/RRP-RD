close all;

%% signal definition
L = 4096;
t = (0:L-1)'/L;

phi1 = 256*t+2700*(t.^2)/2;
phi2 = 768*t+3000*(t.^2)/2;
IF1 = 256+2700*t;
IF2 = 768+3000*t;
s1 = exp(2*1i*pi*phi1);
s2 = exp(2*1i*pi*phi2);
modes_LC = transpose([s1, s2]);
IFs_LC = transpose([IF1, IF2]);
sigma_LC = 0.0188;

Nr = 2;
Nfft = 512;
cas = 1;
poly_degree = 5;

%% SNR (RD + MR) computation
clwin = 10;
SNR_IN = -10:2:4;
NRep = 30;
% SNR_IN = -10:1:-8;
% NRep = 1;

[SNR_NEW_LCR, SNR_NEW_MR, SNR_S_MR, SNR_MB_MR, SNR_IF_NEW, SNR_IF_C_RD, SNR_IF_MB_RD] =...
    test_RD_MR(modes_LC, IFs_LC, clwin, sigma_LC, Nfft, poly_degree, SNR_IN, NRep);

plot_SNR_modes(SNR_IN, SNR_NEW_LCR, SNR_NEW_MR, SNR_S_MR, SNR_MB_MR);
% savefig('F5_MR_LC');
% saveas(gcf,'F5_MR_LC','epsc');
% close all;

plot_SNR_IFs(SNR_IN, SNR_IF_NEW, SNR_IF_C_RD, SNR_IF_MB_RD);
% savefig('F4_RD_LC');
% saveas(gcf,'F4_RD_LC','epsc');
% close all;
