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
A_cos = 2*pi;
phi_s = L/2*t+B/A_cos*cos(A_cos*t);
s_cos = exp(2*1i*pi*phi_s);
sigma_cos = 0.0142;

phip_s = L/2-B*sin(A_cos*t);
phipp_s = -B*A_cos*cos(A_cos*t);

std_s = 1/(sqrt(2*pi)*sigma_cos)*sqrt(1 + sigma_cos^4*phipp_s.^2);
k_min_cos = round((phip_s - std_s)*Nfft/L) + 1;
k_max_cos = round((phip_s + std_s)*Nfft/L) + 1;


%% HOsc
A_osc = 8*pi;
B = floor(2*L/(2*pi));
phi_s = L/2*t+B/A_osc*cos(A_osc*t);
s_osc = exp(2*1i*pi*phi_s);
sigma_osc = 0.0142;

phip_s = L/2-B*sin(A_osc*t);
phipp_s = -B*A_osc*cos(A_osc*t);

std_s = 1/(sqrt(2*pi)*sigma_osc)*sqrt(1 + sigma_osc^4*phipp_s.^2);
k_min_osc = round((phip_s - std_s)*Nfft/L) + 1;
k_max_osc = round((phip_s + std_s)*Nfft/L) + 1;


%% Loops

SNRs = -10:1:-2;
N_snr = length(SNRs);
N_rep = 30;

%% test max

% [N_hit_mean, N_hit_var] = R1_zone_count_rep(s_LC, L, Nfft, sigma_LC, SNRs, N_rep, k_min_LC, k_max_LC);
% [N_hit_mean_cos, N_hit_var_cos] = R1_zone_count_rep(s_cos, L, Nfft, sigma_cos, SNRs, N_rep, k_min_cos, k_max_cos);
% [N_hit_mean_osc, N_hit_var_osc] = R1_zone_count_rep(s_osc, L, Nfft, sigma_osc, SNRs, N_rep, k_min_osc, k_max_osc);
% 
% save('data_fig3_zone_count.mat', 'N_hit_mean', 'N_hit_var',...
%     'N_hit_mean_cos', 'N_hit_var_cos', 'N_hit_mean_osc', 'N_hit_var_osc');

load('data_fig3_zone_count.mat');

%% figures
c1 = [0, 0, 0];
c2 = [0.6350 0.0780 0.1840];
c3 = [0.9290 0.6940 0.1250];

figure;
hold on;
plot(SNRs, N_hit_mean(1, :)/L, 'k-', 'LineWidth', 2,...
    'MarkerSize', 10, 'Color',c1,...
    'DisplayName', '$\mathcal{P}(2)$, linear chirp');
plot(SNRs, N_hit_mean(2, :)/L, 'k-o', 'LineWidth', 2,...
    'MarkerSize', 10, 'Color',c1,...
    'DisplayName', 'Max, linear chirp');
plot(SNRs, N_hit_mean(3, :)/L, 'k-*', 'LineWidth', 2,...
    'MarkerSize', 10, 'Color',c1,...
    'DisplayName', '$\mathcal{Q}(2)$, linear chirp');

plot(SNRs, N_hit_mean_cos(1, :)/L, 'r--', 'LineWidth', 2,...
    'MarkerSize', 10, 'Color',c2,...
    'DisplayName', '$\mathcal{P}(2)$, cosine');
plot(SNRs, N_hit_mean_cos(2, :)/L, 'r--o', 'LineWidth', 2,...
    'MarkerSize', 10, 'Color',c2,...
    'DisplayName', 'Max, cosine');
plot(SNRs, N_hit_mean_cos(3, :)/L, 'r--*', 'LineWidth', 2,...
    'MarkerSize', 10, 'Color',c2,...
    'DisplayName', '$\mathcal{Q}(2)$, cosine');

plot(SNRs, N_hit_mean_osc(1, :)/L, 'g-.', 'LineWidth', 2,...
    'MarkerSize', 10, 'Color',c3,...
    'DisplayName', '$\mathcal{P}(2)$, modulated cosine');
plot(SNRs, N_hit_mean_osc(2, :)/L, 'g-.o', 'LineWidth', 2,...
    'MarkerSize', 10, 'Color',c3,...
    'DisplayName', 'Max, modulated cosine');
plot(SNRs, N_hit_mean_osc(2, :)/L, 'g-.*', 'LineWidth', 2,...
    'MarkerSize', 10, 'Color',c3,...
    'DisplayName', '$\mathcal{Q}(2)$, modulated cosine');
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
savefig('fig_R2_prop_mean');
saveas(gcf,'fig_R2_prop_mean','epsc');
% close all

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
