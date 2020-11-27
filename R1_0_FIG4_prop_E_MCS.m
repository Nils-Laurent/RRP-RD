close all;

L = 4096;
t = (0:L-1)'/L;

%% STFT
Nfft = 512;
Nr = 1;
cas = 1;

fx = (0:Nfft-1)*L/Nfft;

%% signals

phi1 = 768*t+3000*(t.^2)/2;
phi2 = 256*t+2700*(t.^2)/2;
s1 = exp(2*1i*pi*phi1);
s2 = exp(2*1i*pi*phi2);
s_LC = s1 + s2;
sigma_LC = 0.0188;

phi1 = 1400*t+1350/(2*pi)*cos(2*pi*t + pi/2);
phi2 = 3400*t+550/(2*pi)*cos(2*pi*t);
s1_cos_lim = exp(2*1i*pi*phi1);
s2 = exp(2*1i*pi*phi2);
s_cos_lim = s1_cos_lim + s2;
sigma_cos_lim = 0.0175;

B1 = log(2000);
phi1 = 2000*t+exp(B1*t)/B1;
phi2 = 200*t+2300*(t.^2)/2;
s1 = exp(2*1i*pi*phi1);
s2 = exp(2*1i*pi*phi2);
s_exp_LC = s1 + s2;
sigma_exp_LC = 0.0241;

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

P_max = 4;

%% figures

[Prop_E_mean, Prop_E_var] = R1_prop_E_rep(s_LC, L, Nfft, sigma_LC, SNRs, N_rep, P_max);
figure;
hold on;
plot(SNRs, Prop_E_mean(1, :), '-', 'LineWidth', 2,...
    'DisplayName', '$P = 1$');
plot(SNRs, Prop_E_mean(2, :), '--', 'LineWidth', 2,...
    'DisplayName', '$P = 2$');
plot(SNRs, Prop_E_mean(3, :), '-.', 'LineWidth', 2,...
    'DisplayName', '$P = 3$');
plot(SNRs, Prop_E_mean(4, :), ':', 'LineWidth', 2,...
    'DisplayName', '$P = 4$');
hold off;
xlabel('SNRs', 'interpreter', 'latex');
ylabel('$E_{modes}$', 'interpreter', 'latex');
xAX = get(gca,'XAxis');
set(xAX,'FontSize', 26);
yAX = get(gca,'YAxis');
set(yAX,'FontSize', 26);
set(groot, 'defaultLegendInterpreter', 'latex');
lgd = legend('Location', 'southeast');
lgd.FontSize = 24;
pbaspect([1 1 1]);
set(gcf, 'Position',  [0, 0, 1000, 1000]);
savefig('fig_R1_prop_E_MCS_LC');
saveas(gcf,'fig_R1_prop_E_MCS_LC','epsc');
close all

[Prop_E_mean_cos, Prop_E_var_cos] = R1_prop_E_rep(s_cos_lim, L, Nfft, sigma_cos_lim, SNRs, N_rep, P_max);
figure;
hold on;
plot(SNRs, Prop_E_mean_cos(1, :), '-', 'LineWidth', 2,...
    'MarkerSize', 10,...
    'DisplayName', '$P = 1$');
plot(SNRs, Prop_E_mean_cos(2, :), '--', 'LineWidth', 2,...
    'MarkerSize', 10,...
    'DisplayName', '$P = 2$');
plot(SNRs, Prop_E_mean_cos(3, :), '-.', 'LineWidth', 2,...
    'MarkerSize', 10,...
    'DisplayName', '$P = 3$');
plot(SNRs, Prop_E_mean_cos(4, :), ':', 'LineWidth', 2,...
    'MarkerSize', 10,...
    'DisplayName', '$P = 4$');
hold off;
xlabel('SNRs', 'interpreter', 'latex');
ylabel('$E_{modes}$', 'interpreter', 'latex');
xAX = get(gca,'XAxis');
set(xAX,'FontSize', 26);
yAX = get(gca,'YAxis');
set(yAX,'FontSize', 26);
set(groot, 'defaultLegendInterpreter', 'latex');
lgd = legend('Location', 'southeast');
lgd.FontSize = 24;
pbaspect([1 1 1]);
set(gcf, 'Position',  [0, 0, 1000, 1000]);
savefig('fig_R1_prop_E_MCS_cos');
saveas(gcf,'fig_R1_prop_E_MCS_cos','epsc');
close all

[Prop_E_mean_exp, Prop_E_var_exp] = R1_prop_E_rep(s_exp_LC, L, Nfft, sigma_exp_LC, SNRs, N_rep, P_max);
figure;
hold on;
plot(SNRs, Prop_E_mean_exp(1, :), '-', 'LineWidth', 2,...
    'MarkerSize', 10,...
    'DisplayName', '$P = 1$');
plot(SNRs, Prop_E_mean_exp(2, :), '--', 'LineWidth', 2,...
    'MarkerSize', 10,...
    'DisplayName', '$P = 2$');
plot(SNRs, Prop_E_mean_exp(3, :), '-.', 'LineWidth', 2,...
    'MarkerSize', 10,...
    'DisplayName', '$P = 3$');
plot(SNRs, Prop_E_mean_exp(4, :), ':', 'LineWidth', 2,...
    'MarkerSize', 10,...
    'DisplayName', '$P = 4$');
hold off;
xlabel('SNRs', 'interpreter', 'latex');
ylabel('$E_{modes}$', 'interpreter', 'latex');
xAX = get(gca,'XAxis');
set(xAX,'FontSize', 26);
yAX = get(gca,'YAxis');
set(yAX,'FontSize', 26);
set(groot, 'defaultLegendInterpreter', 'latex');
lgd = legend('Location', 'southeast');
lgd.FontSize = 24;
pbaspect([1 1 1]);
set(gcf, 'Position',  [0, 0, 1000, 1000]);
savefig('fig_R1_prop_E_MCS_exp');
saveas(gcf,'fig_R1_prop_E_MCS_exp','epsc');
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
