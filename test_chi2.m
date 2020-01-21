close all

%%
L = 4096;
sigma = 1/sqrt(L);
sigma_noise = pi;
noise = sigma_noise*(randn(L,1)+1i*randn(L,1));

Nfft = 512;
cas = 1;

[g, Lg] = create_gaussian_window(L, Nfft, sigma);
TFR_noise = tfrstft(noise, Nfft, cas, g, Lg);

gamma = median(abs(real(TFR_noise(:))))/0.6745;
C1p = 9.2103;
C10p = 4.6052;
TH = gamma*sqrt(C10p);

N_TH = sum(abs(TFR_noise(:)) > TH);
%truth = sigma_noise*norm(g, 2);

fprintf("ratio = %d\n", 100*N_TH/(L*Nfft));

figure;
imagesc((0:L-1)/L, (L/Nfft)*(1:Nfft), abs(TFR_noise));
set(gca,'ydir','normal');
title("noise");
axis square