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

SNR = inf;
noise = randn(L,1) + 1i*randn(L,1);
s_noise = sigmerge(s_in, noise, SNR);
smooth_p = 1 - 10^(-10);
N_Y = 128;
[IF_FSST, m_FSST, FSST4, IF_new, m_simple, m_LCR, STFT_LCR, STFT] =...
    R1_GW_RRP(s_noise, Fs, Nfft, sigma_L, sigma_Fs, smooth_p, N_Y);

[~, Lh] = create_gaussian_window(Fs, Nfft, sigma_Fs);
% X_cmp = (Lh+1):(L-round(Lh/2));
% X_cmp = Lh:(L-round(Lh/4));
X_cmp = 490:(L-round(Lh/4));
SNR_FSST = snr(sw(X_cmp), m_FSST(X_cmp) - sw(X_cmp));
SNR_S = snr(sw(X_cmp), m_simple(X_cmp) - sw(X_cmp));
SNR_LCR = snr(sw(X_cmp), m_LCR(X_cmp) - sw(X_cmp));
[SNR_FSST, SNR_S, SNR_LCR]
% return;


%% figs
set(groot, 'defaultLegendInterpreter', 'latex');
t_fig = X_cmp/Fs;
X_fig = X_cmp;
STFT2 = STFT(1:N_Y, :);
% N_Y = 128;

%% fig FSST4
figure;
imagesc(t_fig, (0:N_Y-1)*Fs/Nfft, abs(FSST4(1:N_Y, X_fig)));
set(gca,'ydir','normal');
colormap(flipud(gray));
axis square
hold on;
plot(t_fig, IF_FSST(X_fig), 'r', 'Linewidth', 2,...
    'DisplayName', 'FSST4');
hold off;
% legend('FSST4 ridge extraction','Location','northwest');
lgd = legend('Location', 'northwest');
lgd.FontSize = 24;

xlabel('time', 'interpreter', 'latex');
ylabel('frequency', 'interpreter', 'latex');
xAX = get(gca,'XAxis');
set(xAX,'FontSize', 26);
yAX = get(gca,'YAxis');
set(yAX,'FontSize', 26);
pbaspect([1 1 1]);
set(gcf, 'Position',  [0, 0, 1000, 1000])
savefig('fig_R1_GW_FSST');
saveas(gcf,'fig_R1_GW_FSST','epsc');
% close all;

%% fig STFT RD
figure;
imagesc(t_fig, (0:N_Y-1)*Fs/Nfft, abs(STFT2(:, X_fig)));
set(gca,'ydir','normal');
colormap(flipud(gray));
axis square;
hold on;
plot(t_fig, IF_FSST(1, X_fig), 'r--',...
    'Linewidth', 2, 'DisplayName', 'FSST4');
plot(t_fig, IF_new(1, X_fig), 'g',...
    'Linewidth', 2, 'DisplayName', 'RRP-RD, $\lambda = 10^{-10}$');
hold off;
lgd = legend('Location', 'northwest');
lgd.FontSize = 24;

xlabel('time', 'interpreter', 'latex');
ylabel('frequency', 'interpreter', 'latex');
xAX = get(gca,'XAxis');
set(xAX,'FontSize', 26);
yAX = get(gca,'YAxis');
set(yAX,'FontSize', 26);
pbaspect([1 1 1]);
set(gcf, 'Position',  [0, 0, 1000, 1000])
savefig('fig_R1_GW_STFT');
saveas(gcf,'fig_R1_GW_STFT','epsc');
% close all;
% return;

%% fig LCR TFR
figure;
imagesc(t_fig, (0:N_Y-1)*Fs/Nfft, abs(STFT_LCR(:, X_fig)));
set(gca,'ydir','normal');
colormap(flipud(gray));
axis square
xlabel('time', 'interpreter', 'latex');
ylabel('frequency', 'interpreter', 'latex');
xAX = get(gca,'XAxis');
set(xAX,'FontSize', 26);
yAX = get(gca,'YAxis');
set(yAX,'FontSize', 26);
pbaspect([1 1 1]);
set(gcf, 'Position',  [0, 0, 1000, 1000])
savefig('fig_R1_GW_LCR_STFT');
saveas(gcf,'fig_R1_GW_LCR_STFT','epsc');
% close all;

%% reconstructed signals
figure;
hold on;
plot(t_fig, sw(X_fig), 'LineWidth', 2,...
    'DisplayName', 'Numerical relativity');
plot(t_fig, m_LCR(X_fig), '-.', 'LineWidth', 2,...
    'DisplayName', 'RRP-MR-LCR, $\lambda = 10^{-10}$');
hold off;
xlabel('time', 'interpreter', 'latex');
ylabel('amplitude', 'interpreter', 'latex');
% xlim([TCs_vec(3, 1)*T, TCs_vec(3, end)*T]);
xlim([t_fig(1), t_fig(end)]);
lgd = legend;
lgd.FontSize = 24;
xAX = get(gca,'XAxis');
set(xAX,'FontSize', 26);
yAX = get(gca,'YAxis');
set(yAX,'FontSize', 26);
pbaspect([1 1 1]);
set(gcf, 'Position',  [0, 0, 1000, 1000])
savefig('fig_R1_GW_MR');
saveas(gcf,'fig_R1_GW_MR','epsc');
% close all;