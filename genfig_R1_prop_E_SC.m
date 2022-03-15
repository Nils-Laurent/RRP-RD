close all;

L = 4096;
t = (0:L-1)'/L;

%% STFT
Nfft = 512;
Nr = 1;
cas = 1;

fx = (0:Nfft-1)*L/Nfft;

%% signals

B = 3*L/4;
phi_LC = L/8*t+B*(t.^2)/2;
s_LC = exp(2*1i*pi*phi_LC);
sigma_LC = 1/sqrt(B);

B = floor(2*L/(2*pi));
phi_cos = L/2*t+B/(2*pi)*cos(2*pi*t);
s_cos = exp(2*1i*pi*phi_cos);
sigma_cos = 0.0142;

B = log(4096 - 510);
phi_exp = 500*t+exp(B*t)/B;
s_exp = exp(2*1i*pi*phi_exp);
sigma_exp = 0.028;

%% Loops

SNRs = -10:1:0;
% N_snr = length(SNRs);
N_rep = 10;

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

P_max = 3;

%% figures
[Prop_E_mean, Prop_E_var] = R1_prop_E_rep(s_LC, L, Nfft, sigma_LC, SNRs, N_rep, P_max);
figure;
hold on;
plot(SNRs, Prop_E_mean(1, :), '-', 'LineWidth', 2,...
    'DisplayName', '$P = 1$');
plot(SNRs, Prop_E_mean(2, :), ':', 'LineWidth', 2,...
    'DisplayName', '$P = 2$');
plot(SNRs, Prop_E_mean(3, :), '--', 'LineWidth', 2,...
    'DisplayName', '$P = 3$');
hold off;
xlabel('SNRs', 'interpreter', 'latex');
ylabel('Proportion of energy', 'interpreter', 'latex');
xAX = get(gca,'XAxis');
set(xAX,'FontSize', 26);
yAX = get(gca,'YAxis');
set(yAX,'FontSize', 26);
set(groot, 'defaultLegendInterpreter', 'latex');
lgd = legend('Location', 'southeast');
lgd.FontSize = 24;
pbaspect([1 1 1]);
set(gcf, 'Position',  [0, 0, 1000, 1000]);
savefig('fig_R1_prop_E_SC_LC');
saveas(gcf,'fig_R1_prop_E_SC_LC','epsc');
close all

[Prop_E_mean_cos, Prop_E_var_cos] = R1_prop_E_rep(s_cos, L, Nfft, sigma_cos, SNRs, N_rep, P_max);
figure;
hold on;
plot(SNRs, Prop_E_mean_cos(1, :), '-o', 'LineWidth', 2,...
    'MarkerSize', 10,...
    'DisplayName', '$P = 1$');
plot(SNRs, Prop_E_mean_cos(2, :), ':o', 'LineWidth', 2,...
    'MarkerSize', 10,...
    'DisplayName', '$P = 2$');
plot(SNRs, Prop_E_mean_cos(3, :), '--o', 'LineWidth', 2,...
    'MarkerSize', 10,...
    'DisplayName', '$P = 3$');
hold off;
xlabel('SNRs', 'interpreter', 'latex');
ylabel('Proportion of energy', 'interpreter', 'latex');
xAX = get(gca,'XAxis');
set(xAX,'FontSize', 26);
yAX = get(gca,'YAxis');
set(yAX,'FontSize', 26);
set(groot, 'defaultLegendInterpreter', 'latex');
lgd = legend('Location', 'southeast');
lgd.FontSize = 24;
pbaspect([1 1 1]);
set(gcf, 'Position',  [0, 0, 1000, 1000]);
savefig('fig_R1_prop_E_SC_cos');
saveas(gcf,'fig_R1_prop_E_SC_cos','epsc');
close all

[Prop_E_mean_exp, Prop_E_var_exp] = R1_prop_E_rep(s_exp, L, Nfft, sigma_exp, SNRs, N_rep, P_max);
figure;
hold on;
plot(SNRs, Prop_E_mean_exp(1, :), '-o', 'LineWidth', 2,...
    'MarkerSize', 10,...
    'DisplayName', '$P = 1$');
plot(SNRs, Prop_E_mean_exp(2, :), ':o', 'LineWidth', 2,...
    'MarkerSize', 10,...
    'DisplayName', '$P = 2$');
plot(SNRs, Prop_E_mean_exp(3, :), '--o', 'LineWidth', 2,...
    'MarkerSize', 10,...
    'DisplayName', '$P = 3$');
hold off;
xlabel('SNRs', 'interpreter', 'latex');
ylabel('Proportion of energy', 'interpreter', 'latex');
xAX = get(gca,'XAxis');
set(xAX,'FontSize', 26);
yAX = get(gca,'YAxis');
set(yAX,'FontSize', 26);
set(groot, 'defaultLegendInterpreter', 'latex');
lgd = legend('Location', 'southeast');
lgd.FontSize = 24;
pbaspect([1 1 1]);
set(gcf, 'Position',  [0, 0, 1000, 1000]);
savefig('fig_R1_prop_E_SC_exp');
saveas(gcf,'fig_R1_prop_E_SC_exp','epsc');
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
