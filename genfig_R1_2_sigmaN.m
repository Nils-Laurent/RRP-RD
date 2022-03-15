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
%% fig STFT
% [g, Lh] = create_gaussian_window(L, Nfft, sigma_LC);
% [STFT] = tfrstft(s_LC, Nfft, 1, g, Lh);
% 
% figure;
% imagesc(t, fx, abs(STFT));
% set(gca,'ydir','normal');
% colormap(flipud(gray));
% axis square;
% % xlabel('time', 'interpreter', 'latex');
% % ylabel('frequency', 'interpreter', 'latex');
% yticks([]);
% xticks([]);
% xAX = get(gca,'XAxis');
% set(xAX,'FontSize', 26);
% yAX = get(gca,'YAxis');
% set(yAX,'FontSize', 26);
% pbaspect([1 0.3 1]);
% set(gcf, 'Position',  [0, 0, 1000, 300])
% savefig('fig_R2_sigma_LC_stft');
% saveas(gcf,'fig_R2_sigma_LC_stft','epsc');
% close all
% 
% 
% [g, Lh] = create_gaussian_window(L, Nfft, sigma_cos);
% [STFT_cos] = tfrstft(s_cos, Nfft, 1, g, Lh);
% 
% figure;
% imagesc(t, fx, abs(STFT_cos));
% set(gca,'ydir','normal');
% colormap(flipud(gray));
% axis square;
% % xlabel('time', 'interpreter', 'latex');
% % ylabel('frequency', 'interpreter', 'latex');
% yticks([]);
% xticks([]);
% xAX = get(gca,'XAxis');
% set(xAX,'FontSize', 26);
% yAX = get(gca,'YAxis');
% set(yAX,'FontSize', 26);
% pbaspect([1 0.3 1]);
% set(gcf, 'Position',  [0, 0, 1000, 300])
% savefig('fig_R2_sigma_cos_stft');
% saveas(gcf,'fig_R2_sigma_cos_stft','epsc');
% close all
% 
% 
% [g, Lh] = create_gaussian_window(L, Nfft, sigma_osc);
% [STFT_osc] = tfrstft(s_osc, Nfft, 1, g, Lh);
% 
% figure;
% imagesc(t, fx, abs(STFT_osc));
% set(gca,'ydir','normal');
% colormap(flipud(gray));
% axis square;
% % xlabel('time', 'interpreter', 'latex');
% % ylabel('frequency', 'interpreter', 'latex');
% yticks([]);
% xticks([]);
% xAX = get(gca,'XAxis');
% set(xAX,'FontSize', 26);
% yAX = get(gca,'YAxis');
% set(yAX,'FontSize', 26);
% pbaspect([1 0.3 1]);
% set(gcf, 'Position',  [0, 0, 1000, 300])
% savefig('fig_R2_sigma_osc_stft');
% saveas(gcf,'fig_R2_sigma_osc_stft','epsc');
% close all

% return;
%% Loops

% SNRs = [-10, -8, -5, -2, 0];
SNR = -10;
N_SNRs = length(SNR);
Y_th = 1.0:0.1:3;
N_th = length(Y_th);
N_rep = 100;

% %% count LC
% [c_global_LC, counts_LC] = R1_global_count_rep(s_LC, L, Nfft, sigma_LC, SNR, Y_th, N_rep, k_min_LC, k_max_LC);
% 
% %% count cos
% [c_global_cos, counts_cos] = R1_global_count_rep(s_cos, L, Nfft, sigma_cos, SNR, Y_th, N_rep, k_min_cos, k_max_cos);
% 
% %% count osc
% [c_global_osc, counts_osc] = R1_global_count_rep(s_osc, L, Nfft, sigma_osc, SNR, Y_th, N_rep, k_min_osc, k_max_osc);
% 
% 
% save('data_fig2_sigmaN.mat', 'c_global_LC', 'counts_LC',...
%     'c_global_cos', 'counts_cos', 'c_global_osc', 'counts_osc');

load('data_fig2_sigmaN.mat');

N_hit_mean = mean(counts_LC, 3);
N_hit_std = std(counts_LC, 0, 3);
N_hit_global_mean = mean(c_global_LC, 3);
N_hit_global_std = std(c_global_LC, 0, 3);

N_hit_mean_cos = mean(counts_cos, 3);
N_hit_std_cos = std(counts_cos, 0, 3);
N_hit_global_mean_cos = mean(c_global_cos, 3);
N_hit_global_std_cos = std(c_global_cos, 0, 3);

N_hit_mean_osc = mean(counts_osc, 3);
N_hit_std_osc = std(counts_osc, 0, 3);
N_hit_global_mean_osc = mean(c_global_osc, 3);
N_hit_global_std_osc = std(c_global_osc, 0, 3);

%% figs LC
figure;
hold on;
plot(Y_th, N_hit_mean(1, :)/L, 'LineWidth', 2,...
    'DisplayName', "$\mathcal{P}(\beta)$, linear chirp");
plot(Y_th, N_hit_mean_cos(1, :)/L, '--', 'LineWidth', 2,...
    'DisplayName', "$\mathcal{P}(\beta)$, cosine");
plot(Y_th, N_hit_mean_osc(1, :)/L, '-.', 'LineWidth', 2,...
    'DisplayName', "$\mathcal{P}(\beta)$, modulated cosine");
plot(Y_th, N_hit_global_mean(1, :)/L, '-o', 'LineWidth', 2,...
    'DisplayName', "Max, linear chirp");
plot(Y_th, N_hit_global_mean_cos(1, :)/L, '--o', 'LineWidth', 2,...
    'DisplayName', "Max, cosine");
plot(Y_th, N_hit_global_mean_osc(1, :)/L, '-.o', 'LineWidth', 2,...
    'DisplayName', "Max, modulated cosine");
hold off;
xlim([Y_th(1), Y_th(end)]);
xlabel('$\beta$', 'interpreter', 'latex');
ylabel('Proportion of detection', 'interpreter', 'latex');
xAX = get(gca,'XAxis');
set(xAX,'FontSize', 24);
yAX = get(gca,'YAxis');
set(yAX,'FontSize', 24);
set(groot, 'defaultLegendInterpreter', 'latex');
lgd = legend('Location', 'southwest');
lgd.FontSize = 24;
pbaspect([1 1 1]);
set(gcf, 'Position',  [0, 0, 1000, 1000])
savefig('fig_R2_sigma_mean');
saveas(gcf,'fig_R2_sigma_mean','epsc');
% close all

% figure;
% hold on;
% plot(Y_th, N_hit_std(1, :), 'LineWidth', 2,...
%     'DisplayName', "linear chirp");
% plot(Y_th, N_hit_std_cos(1, :), '--', 'LineWidth', 2,...
%     'DisplayName', "cosine");
% plot(Y_th, N_hit_global_std(1, :), '-o', 'LineWidth', 2,...
%     'DisplayName', "Max, linear chirp");
% plot(Y_th, N_hit_global_std_cos(1, :), '--o', 'LineWidth', 2,...
%     'DisplayName', "Max, cosine");
% hold off;
% xlim([Y_th(1), Y_th(end)]);
% xlabel('$\beta$', 'interpreter', 'latex');
% ylabel('std', 'interpreter', 'latex');
% xAX = get(gca,'XAxis');
% set(xAX,'FontSize', 24);
% yAX = get(gca,'YAxis');
% set(yAX,'FontSize', 24);
% set(groot, 'defaultLegendInterpreter', 'latex');
% lgd = legend('Location', 'northwest');
% lgd.FontSize = 24;
% pbaspect([1 1 1]);
% set(gcf, 'Position',  [0, 0, 1000, 1000])
% savefig('fig_R1_sigma_std');
% saveas(gcf,'fig_R1_sigma_std','epsc');
% close all
