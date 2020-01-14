close all

%% signal definition
N = 4096;
t = (0:N-1)/N;

s = 2.*exp(2*pi*1i*(1000*t+30*cos(6*pi*t)));
s = s(:);

%% apply noise

var = 1;
WGN = var*randn(N,1)+var*1i*randn(N,1);
s_clean = s;
s = sigmerge(s, WGN, -10);

%% STFT
Nfft = 512;
downsampling = 1;
N_down = N/downsampling;
shift = 0;
cas = 1;

WLengths = [0.04, 0.055, 0.15, 0.35, 0.5];
NWin = length(WLengths);
TFRs = zeros(Nfft, N_down, NWin);
for k = 1:NWin
    [h, Lh] = tftb_window2('gauss', Nfft, WLengths(k));
    [TFRs(:, :, k), norm2h] = tfrstft_three_case_down(s, Nfft, cas, h, Lh, downsampling, shift);
    disp(["k STFT ", k]);
end

%% hard thresholding
disp("hard thresholding");
TFRs_d = TFRs;
for k = 1:NWin
    y2 = real(TFRs_d(:, :, k));
    th = median(abs(y2(:)))/0.6745;
    th_max = 0.7*max(abs(TFRs_d(:, :, k)));
    
    TFRk = TFRs_d(:, :, k);
    TFRk(abs(TFRk) < th_max) = 0;
    TFRs_d(:, :, k) = TFRk;
end

%% ridge extraction
disp("ridge extraction");
%[Cs_d] = exridge(tfr_d,0,0,50);
 CSk = zeros(1, N_down, NWin);
for k = 1:NWin
    [Cs] = exridge(TFRs_d(:, :, k),0,0,50);
    Cs = Cs';
    CSk(:, :, k) = Cs;
end

%% display result
disp("display result");
% figure;
% for k = 1:NWin
%     subplot(2, 5, k);
%     imagesc((0:N_down-1)/N_down, 1:Nfft, abs(TFRs_d(:, :, k)));
%     set(gca,'ydir','normal');
%     title(["winL = ", WLengths(k)]);
%     axis square
%     
%     subplot(2, 5, k + NWin);
%     imagesc((0:N_down-1)/N_down, 1:Nfft, abs(TFRs_d(:, :, k)));
%     set(gca,'ydir','normal');
%     title(["winL = ", WLengths(k)]);
%     hold on;
%     plot((0:N_down-1)/N_down, CSk(1, :, k)-1,'r','linewidth',1);
%     hold off;
%     axis square
% end

%figure;
for k = 1:NWin
    %subplot(2, ceil(NWin/2), k);
    figure;
    imagesc((0:N_down-1)/N_down, 1:Nfft, abs(TFRs_d(:, :, k)));
    set(gca,'ydir','normal');
    title(["winL = ", WLengths(k)]);
    axis square
end

figure;
%subplot(1, 2, 1);
TFRsum = TFRs_d(:, :, 1);
for k = 2:NWin
    TFRsum = TFRsum + TFRs_d(:, :, k);
end

imagesc((0:N_down-1)/N_down, 1:Nfft, abs(TFRsum));
set(gca,'ydir','normal');
title("SUM");
axis square