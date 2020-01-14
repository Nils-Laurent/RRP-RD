% Analysis case : Ridge detection with high noise level

close all

%% signal definition
L = 4096;
t = (0:L-1)'/L;

A = 2048;
B = 3/4*L/2;
phi = A*t+B*cos(4*pi*t)/(4*pi);
s_clean = exp(2*1i*pi*phi);

%% apply noise

var = 1;
WGN = var*randn(L,1)+var*1i*randn(L,1);
s_noise = sigmerge(s_clean, WGN, -10);

%% STFT
Nfft = 512;
%sigma = 1/sqrt(4*pi*B);
sigma = 1/((7/3)^(1/4)*sqrt(4*pi*B));
cas = 1;

[g, Lg] = create_gaussian_window(L, Nfft, sigma);
[TFR_clean] = tfrstft(s_clean, Nfft, cas, g, Lg);
[TFR_noise] = tfrstft(s_noise, Nfft, cas, g, Lg);

[maxv, arg] = max(abs(TFR_noise(:)));

clwin = Nfft;
[~, TFR_denoised, Cs] = denoise_LCD(s_noise, 1, clwin, sigma, Nfft, 0, 0);

%% display

figure;
subplot(1, 3, 1);
imagesc((0:L-1)/L, (L/Nfft)*(1:Nfft), abs(TFR_clean));
set(gca,'ydir','normal');
title("TFR clean");
axis square

% figure;
subplot(1, 3, 2);
imagesc((0:L-1)/L, (L/Nfft)*(1:Nfft), abs(TFR_noise));
set(gca,'ydir','normal');
title("TFR noise + RD");
axis square

hold on;
argX = floor(arg/Nfft);
argY = mod(arg, Nfft);
plot((argX - 1)/L, (L/Nfft)*(argY), '-o', 'MarkerSize', 15, 'MarkerEdgeColor','red');
plot((0:L-1)/L, (L/Nfft)*(Cs-1),'r','linewidth',1);
hold off;

% figure;
subplot(1, 3, 3);
imagesc((0:L-1)/L, (L/Nfft)*(1:Nfft), abs(TFR_denoised));
set(gca,'ydir','normal');
title("TFR denoised");
axis square

hold on;
argX = floor(arg/Nfft);
argY = mod(arg, Nfft);
plot((argX - 1)/L, (L/Nfft)*(argY), '-o', 'MarkerSize', 15, 'MarkerEdgeColor','red');
hold off;
