close all;

load GW-observed-Hanford.txt
t1=GW_observed_Hanford(:,1);
s1=GW_observed_Hanford(:,2);
% figure;
% hold on;
% plot(t1(1:end-1));
% plot(t1(2:end));
% plot(tw, '--');
% hold off;
% pause;
t1 = t1(2:end);
s1 = s1(2:end);

load GW-waveform-H.txt
tw=GW_waveform_H(:,1);
sw=GW_waveform_H(:,2)';

t = t1;
s_in = s1(:);
L = length(s_in);
s_in = hilbert(s_in);
T = 0.21;
Fs = round(L/T);
sigma_L = 0.05;
sigma_Fs = sigma_L*L/Fs;
Nfft = 4096;
Nr = 1;


SNRs_IN = -18:1:0;
N_SNR = length(SNRs_IN);
smooth_p = 1 - 10^(-10);
N_Y = 128;

[~, Lh] = create_gaussian_window(Fs, Nfft, sigma_Fs);
% X_cmp = (Lh+1):(L-round(Lh/2));
X_cmp = 490:(L-round(Lh/4));
% n_rd = 2965;
% X1 = 1425:n_rd;
% X2 = n_rd:(L-round(Lh/4));

RES_SNR = zeros(3, N_SNR);
N_rep = 40;
% load('snr_R1_GW.mat');
for n=1:N_SNR
%     break;
    for m=1:N_rep
        fprintf("snr = %d, rep = %d\n", SNRs_IN(n), m);
        SNR = SNRs_IN(n);
        noise = randn(L,1) + 1i*randn(L,1);
        s_noise = sigmerge(s_in, noise, SNR);
        [IF_FSST, m_FSST, FSST4, IF_new, m_simple, m_LCR, STFT_LCR, STFT] =...
            R1_GW_RRP(s_noise, Fs, Nfft, sigma_L, sigma_Fs, smooth_p, N_Y);

        X_set = X_cmp;
        c_FSST = snr(sw(X_set), m_FSST(X_set) - sw(X_set));
        c_simple = snr(sw(X_set), m_simple(X_set) - sw(X_set));
        c_LCR = snr(sw(X_set), m_LCR(X_set) - sw(X_set));
        RES_SNR(1, n) = RES_SNR(1, n) + c_FSST/N_rep;
        RES_SNR(2, n) = RES_SNR(2, n) + c_simple/N_rep;
        RES_SNR(3, n) = RES_SNR(3, n) + c_LCR/N_rep;
    end
    
    %% info
%     fprintf("--- SNR = %d ---\n", SNR);
%     RES_SNR(:, n)
%     
%     t_fig = X_cmp/Fs;
%     X_fig = X_cmp;
%     figure;
%     imagesc(t_fig, (0:N_Y-1)*Fs/Nfft, abs(STFT(1:N_Y, X_fig)));
%     set(gca,'ydir','normal');
%     colormap(flipud(gray));
%     axis square;
%     hold on;
%     plot(t_fig, IF_new(1, X_fig));
%     plot(t_fig, IF_FSST(1, X_fig), '--');
%     hold off;
%     pause;
end
save('snr_R1_GW.mat', 'RES_SNR');


%% figs
set(groot, 'defaultLegendInterpreter', 'latex');

B1 = [0 0.4470 0.7410]/2;
B3 = (1 + [0.3010 0.7450 0.9330])/2;
B2 = (B3 + B1)/2;

R_Cl = [0.6350 0.0780 0.1840];
P_Cl = [0.4940 0.1840 0.5560];
G_FSST = [0.4660 0.6740 0.1880];
R2 = [0.8500 0.3250 0.0980];

figure;
hold on;
plot(SNRs_IN, RES_SNR(2, :), '--',...
    'Linewidth', 2, 'Markersize', 10,...
    'DisplayName', 'RRP-MR, $\lambda = 10^{-10}$',...
    'Color', B1);
plot(SNRs_IN, RES_SNR(1, :), '-o',...
    'Linewidth', 2, 'Markersize', 10,...
    'DisplayName', 'FSST4-MR',...
    'Color', B2);
%     'Color', R_MB);
plot(SNRs_IN, RES_SNR(3, :), '-*',...
    'Linewidth', 2, 'Markersize', 10,...
    'DisplayName', 'RRP-LCR-MR, $\lambda = 10^{-10}$',...
    'Color', R2);
hold off;
xlim([SNRs_IN(1), SNRs_IN(end)]);
lgd = legend('Location', 'southeast');
lgd.FontSize = 24;
xlabel('SNR in', 'interpreter', 'latex');
ylabel('SNR out', 'interpreter', 'latex');
xAX = get(gca,'XAxis');
set(xAX,'FontSize', 26);
yAX = get(gca,'YAxis');
set(yAX,'FontSize', 26);
pbaspect([1 1 1]);
set(gcf, 'Position',  [0, 0, 1000, 1000])
savefig('fig_R1_GW_SNR');
saveas(gcf, 'fig_R1_GW_SNR', 'epsc');