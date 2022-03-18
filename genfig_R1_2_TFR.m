close all;

conf_genfig;

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
[g, ~] = gauss_win(L, sigma_LC);
[STFT, ~] = stft(s_LC, Nfft, g);

figure;
imagesc(t, fx, abs(STFT));
set(gca,'ydir','normal');
colormap(flipud(gray));
axis square;
% xlabel('time', 'interpreter', 'latex');
% ylabel('frequency', 'interpreter', 'latex');
yticks([]);
xticks([]);
xAX = get(gca,'XAxis');
set(xAX,'FontSize', 26);
yAX = get(gca,'YAxis');
set(yAX,'FontSize', 26);
pbaspect([1 0.3 1]);
set(gcf, 'Position',  [0, 0, 1000, 300])
savefig('fig_R2_sigma_LC_stft');
saveas(gcf,'fig_R2_sigma_LC_stft','epsc');
close all


[g, ~] = gauss_win(L, sigma_cos);
[STFT_cos, ~] = stft(s_cos, Nfft, g);

figure;
imagesc(t, fx, abs(STFT_cos));
set(gca,'ydir','normal');
colormap(flipud(gray));
axis square;
% xlabel('time', 'interpreter', 'latex');
% ylabel('frequency', 'interpreter', 'latex');
yticks([]);
xticks([]);
xAX = get(gca,'XAxis');
set(xAX,'FontSize', 26);
yAX = get(gca,'YAxis');
set(yAX,'FontSize', 26);
pbaspect([1 0.3 1]);
set(gcf, 'Position',  [0, 0, 1000, 300])
savefig('fig_R2_sigma_cos_stft');
saveas(gcf,'fig_R2_sigma_cos_stft','epsc');
close all


[g, ~] = gauss_win(L, sigma_osc);
[STFT_osc, ~] = stft(s_osc, Nfft, g);

figure;
imagesc(t, fx, abs(STFT_osc));
set(gca,'ydir','normal');
colormap(flipud(gray));
axis square;
% xlabel('time', 'interpreter', 'latex');
% ylabel('frequency', 'interpreter', 'latex');
yticks([]);
xticks([]);
xAX = get(gca,'XAxis');
set(xAX,'FontSize', 26);
yAX = get(gca,'YAxis');
set(yAX,'FontSize', 26);
pbaspect([1 0.3 1]);
set(gcf, 'Position',  [0, 0, 1000, 300])
savefig('fig_R2_sigma_osc_stft');
saveas(gcf,'fig_R2_sigma_osc_stft','epsc');
close all

