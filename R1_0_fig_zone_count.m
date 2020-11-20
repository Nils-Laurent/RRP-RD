close all;

L = 4096;
t = (0:L-1)'/L;

%% STFT
Nfft = 512;
Nr = 1;
cas = 1;

fx = (0:Nfft-1)*L/Nfft;

%% linear chirp signal
A = L/32;
B = L - 2*A;

phi_s = A*t+B/2*(t.^2);
s_LC = exp(2*1i*pi*phi_s);
sigma_LC = 1/sqrt(B);

phip_s = A + B*t;
phipp_s = B*ones(L, 1);

std_s = 1/(sqrt(2*pi)*sigma_LC)*sqrt(1 + sigma_LC^4*phipp_s.^2);
k_min_LC = round((phip_s - std_s)*Nfft/L) + 1;
k_max_LC = round((phip_s + std_s)*Nfft/L) + 1;


%% cosine signal
B = floor(2*L/(2*pi));
phi_s = L/2*t+B/(2*pi)*cos(2*pi*t);
s_cos = exp(2*1i*pi*phi_s);
sigma_cos = 0.0142;

phip_s = L/2-B*sin(2*pi*t);
phipp_s = -2*pi*cos(2*pi*t);

std_s = 1/(sqrt(2*pi)*sigma_cos)*sqrt(1 + sigma_cos^4*phipp_s.^2);
k_min_cos = round((phip_s - std_s)*Nfft/L) + 1;
k_max_cos = round((phip_s + std_s)*Nfft/L) + 1;

%% Loops

SNRs = -10:1:-2;
% N_snr = length(SNRs);
N_rep = 30;

% y = 2;

%% test max
% SNR = -10;
% noise = randn(L,1) + 1i*randn(L,1);
% % load('noise_frag.mat');
% s_noise = sigmerge(s_in, noise, SNR);
% [g, Lh] = create_gaussian_window(L, Nfft, sigma_s);
% [STFT, omega, ~, QM, ~, tau] = FM_operators(s_noise, Nfft, g, Lh, sigma_s);
% R1_TF_in_STD(STFT, QM, omega, tau, L, Nfft, y, k_IF_min, k_IF_max)
% return;

[N_hit_mean, N_hit_var] = R1_zone_count_rep(s_LC, L, Nfft, sigma_LC, SNRs, N_rep, k_min_LC, k_max_LC);
[N_hit_mean_cos, N_hit_var_cos] = R1_zone_count_rep(s_cos, L, Nfft, sigma_cos, SNRs, N_rep, k_min_cos, k_max_cos);

%% figures
figure;
hold on;
plot(SNRs, N_hit_mean(1, :)/L, '-', 'LineWidth', 2,...
    'DisplayName', '$\mathcal{P}(2)$, linear chirp');
plot(SNRs, N_hit_mean(2, :)/L, ':', 'LineWidth', 2,...
    'DisplayName', 'Max, linear chirp');
plot(SNRs, N_hit_mean(3, :)/L, '--', 'LineWidth', 2,...
    'DisplayName', '$\mathcal{Q}(2)$, linear chirp');
plot(SNRs, N_hit_mean_cos(1, :)/L, '-o', 'LineWidth', 2,...
    'MarkerSize', 10,...
    'DisplayName', '$\mathcal{P}(2)$, cosine');
plot(SNRs, N_hit_mean_cos(2, :)/L, ':o', 'LineWidth', 2,...
    'MarkerSize', 10,...
    'DisplayName', 'Max, cosine');
plot(SNRs, N_hit_mean_cos(3, :)/L, '--o', 'LineWidth', 2,...
    'MarkerSize', 10,...
    'DisplayName', '$\mathcal{Q}(2)$, cosine');
hold off;
xlabel('SNRs', 'interpreter', 'latex');
ylabel('Proportion of detection', 'interpreter', 'latex');
xAX = get(gca,'XAxis');
set(xAX,'FontSize', 26);
yAX = get(gca,'YAxis');
set(yAX,'FontSize', 26);
set(groot, 'defaultLegendInterpreter', 'latex');
lgd = legend('Location', 'southeast');
lgd.FontSize = 24;
pbaspect([1 1 1]);
set(gcf, 'Position',  [0, 0, 1000, 1000]);
savefig('fig_R1_prop_mean');
saveas(gcf,'fig_R1_prop_mean','epsc');
close all

% figure;
% hold on;
% plot(SNRs, N_hit_mean_cos(1, :), 'LineWidth', 2,...
%     'DisplayName', '$\mathcal{P}(\beta)$');
% plot(SNRs, N_hit_mean_cos(2, :), 'g--o', 'LineWidth', 2,...
%     'DisplayName', 'Max');
% plot(SNRs, N_hit_mean_cos(3, :), 'r--', 'LineWidth', 2,...
%     'DisplayName', 'Zone max');
% hold off;
% xlabel('SNRs', 'interpreter', 'latex');
% ylabel('Proportion of detection', 'interpreter', 'latex');
% xAX = get(gca,'XAxis');
% set(xAX,'FontSize', 26);
% yAX = get(gca,'YAxis');
% set(yAX,'FontSize', 26);
% set(groot, 'defaultLegendInterpreter', 'latex');
% lgd = legend('Location', 'southeast');
% lgd.FontSize = 24;
% pbaspect([1 1 1]);
% set(gcf, 'Position',  [0, 0, 1000, 1000]);
% savefig('fig_R1_prop_mean_cos');
% saveas(gcf,'fig_R1_prop_mean_cos','epsc');
% close all

% figure;
% hold on;
% plot(SNRs, N_hit_var(1, :), 'LineWidth', 2,...
%     'DisplayName', '$\mathcal{P}(\beta)$');
% plot(SNRs, N_hit_var(2, :), 'g--o', 'LineWidth', 2,...
%     'DisplayName', 'Max');
% plot(SNRs, N_hit_var(3, :), 'r--', 'LineWidth', 2,...
%     'DisplayName', 'Zone max');
% hold off;
% xlabel('SNRs', 'interpreter', 'latex');
% ylabel('std', 'interpreter', 'latex');
% xAX = get(gca,'XAxis');
% set(xAX,'FontSize', 26);
% yAX = get(gca,'YAxis');
% set(yAX,'FontSize', 26);
% set(groot, 'defaultLegendInterpreter', 'latex');
% lgd = legend('Location', 'southeast');
% lgd.FontSize = 24;
% pbaspect([1 1 1]);
% set(gcf, 'Position',  [0, 0, 1000, 1000]);
% savefig('fig_R1_prop_var');
% saveas(gcf,'fig_R1_prop_var','epsc');
% close all
