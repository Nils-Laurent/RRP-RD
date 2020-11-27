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
% return;

% [g, Lg] = create_gaussian_window(L, Nfft, sigma_s);
% [STFT, ~, ~, QM] = FM_operators(s, Nfft, g, Lg, sigma_s);

% SNR = -3;
% noise = randn(L,1) + 1i*randn(L,1);
% s_noise = sigmerge(s_in, noise, SNR);
s_noise = s_in;

%% FSST4
% 1: verif meme fenetre

Te1=t1(2)-t1(1);
nv = log2(L);
if mod(nv,1)~=0
    warning('The signal is not a power of two, zero padding to the next power');
    s1_pad = [s1.' zeros(1,2^(floor(log2(L))+1)-L)];
    t1pad = max(t1)+Te1:Te1:max(t1)+(2^(floor(log2(L))+1)-L)*Te1;
    t1_pad = [t1.' t1pad];
    sw_pad = [sw zeros(1,2^(floor(log2(L))+1)-L)];
end
L_pad = length(s1_pad);
s_pad = hilbert(s1_pad);

gamma = 10^(-3);
clwin = 10;
lambda = 0;
% L_pad = L;
% s_pad = s_noise;
ft =1:L_pad/32;bt=1:L_pad;

[STFT,FSST,FSST2,FSST3,FSST4,omega,omega2,omega3,tau2,tau3,phi22p,phi33p,phi44p] =...
    sstn(s_pad,gamma,sigma_L,ft,bt);
[Cs, Es] = exridge_mult_Noise(FSST4, Nr, lambda, clwin);

%% RRP
Fs1 = Fs;
% Fs = L;
sigma_s = sigma_Fs;

[g, Lh] = create_gaussian_window(Fs, Nfft, sigma_s);
% X_cmp = Lh+1:(L-Lh);
X_cmp = 688:(L-Lh);
[STFT, omega, ~, QM, ~, tau] = FM_operators(s_noise, Fs, Nfft, g, Lh, sigma_s);

N_Y = 128;
STFT2 = STFT(1:N_Y, :);
QM2 = QM(1:N_Y, :);
omega2 = omega(1:N_Y, :);
tau2 = tau(1:N_Y, :);

% smooth_vec = [1 - 10^(-3), 1 - 10^(-6), 1 - 10^(-10)];
smooth_vec = 1 - 10^(-6);
Ns = length(smooth_vec);

SNRs = zeros(Ns, 2);
IFs_vec = zeros(Ns, L);

ind = 0;
m_LCR = [];
STFT_LCR = [];
for smooth_p = smooth_vec
    ind = ind + 1;

    [Spl, ~] = R1_RRP_RD(STFT2, QM2, omega2, tau2, Fs, Nfft, Nr, sigma_s, smooth_p);

    [m_simple, m_LCR, IFs_p, STFT_LCR] = R1_MR_and_LCR_spl(STFT2, Spl, g, Lh, sigma_s, Nr, Nfft, Fs);
    IFs_vec(ind, :) = IFs_p(1, :);
    m_simple = real(m_simple);
    m_LCR = real(m_LCR);

    % 2: Retier la demie fenetre
    SNRs(ind, 1) = snr(sw(X_cmp), m_simple(X_cmp) - sw(X_cmp));
    SNRs(ind, 2) = snr(sw(X_cmp), m_LCR(X_cmp) - sw(X_cmp));
    figure;
    hold on;
    plot(X_cmp, sw(X_cmp), 'LineWidth', 2,...
        'DisplayName', 'Numerical relativity');
    plot(X_cmp, m_LCR(X_cmp), '-.', 'LineWidth', 2,...
        'DisplayName', 'RRP-MR-LCR');
    hold off;
    legend;
end
SNRs
IFs_vec = IFs_vec*Fs1/Fs;
Fs = Fs1;
return;


%% figs
set(groot, 'defaultLegendInterpreter', 'latex');
t_fig = X_cmp/Fs;
X_fig = X_cmp;

%% fig FSST4
figure;
imagesc(t_fig, (0:N_Y-1)*Fs/Nfft, abs(FSST4(1:N_Y, X_fig)));
set(gca,'ydir','normal');
colormap(flipud(gray));
axis square
hold on;
plot(t_fig, (Cs(X_fig)-1)*Fs/Nfft, 'g', 'Linewidth', 2);
hold off;
legend('FSST4 ridge extraction','Location','northwest');
% pbaspect([1 1 1]);
% set(gcf, 'Position',  [0, 0, 1000, 1000])

%% fig STFT RD
figure;
imagesc(t_fig, (0:N_Y-1)*Fs/Nfft, abs(STFT2(:, X_fig)));
set(gca,'ydir','normal');
colormap(flipud(gray));
axis square;
hold on;
% plot((0:L-1)/Fs, fnval(Spl(1).spline, (0:L-1)/Fs));
for m=1:Ns
    plot(t_fig, IFs_vec(m, X_fig));
end
hold off;

xlabel('time', 'interpreter', 'latex');
ylabel('frequency', 'interpreter', 'latex');
xAX = get(gca,'XAxis');
set(xAX,'FontSize', 26);
yAX = get(gca,'YAxis');
set(yAX,'FontSize', 26);
% pbaspect([1 1 1]);
% set(gcf, 'Position',  [0, 0, 1000, 1000])
% savefig('F_GW_STFT');
% saveas(gcf,'F_GW_STFT','epsc');
% close all;


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
% pbaspect([1 1 1]);
% set(gcf, 'Position',  [0, 0, 1000, 1000])
% savefig('F_GW_LCR_STFT');
% saveas(gcf,'F_GW_LCR_STFT','epsc');
% close all;

%% reconstructed signals
figure;
hold on;
plot(t_fig, sw(X_fig), 'LineWidth', 2,...
    'DisplayName', 'Numerical relativity');
plot(t_fig, m_LCR(X_fig), '-.', 'LineWidth', 2,...
    'DisplayName', 'RRP-MR-LCR');
hold off;
xlabel('time', 'interpreter', 'latex');
ylabel('dB', 'interpreter', 'latex');
% xlim([TCs_vec(3, 1)*T, TCs_vec(3, end)*T]);
lgd = legend;
lgd.FontSize = 37;
xAX = get(gca,'XAxis');
set(xAX,'FontSize', 26);
yAX = get(gca,'YAxis');
set(yAX,'FontSize', 26);
% pbaspect([1 1 1]);
% set(gcf, 'Position',  [0, 0, 1000, 1000])
% savefig('F_GW_NR');
% saveas(gcf,'F_GW_NR','epsc');
% close all;

%% TFR inter
% figure;
% imagesc((0:L-1)/L*T, (0:127)*L/(Nfft*T), abs(TFR_inter));
% set(gca,'ydir','normal');
% colormap(flipud(gray));
% axis square
% xlabel('time', 'interpreter', 'latex');
% ylabel('frequency', 'interpreter', 'latex');
% hold on;
% plot(TCs_vec(3, :)*T, IFs_vec(3, :)/T,...
%     'LineWidth', 2, 'DisplayName', 'RRP-RD (tol = 3)');
% 
% hold off;
% lgd = legend;
% lgd.FontSize = 37;
% xAX = get(gca,'XAxis');
% set(xAX,'FontSize', 26);
% yAX = get(gca,'YAxis');
% set(yAX,'FontSize', 26);
% pbaspect([1 1 1]);
% set(gcf, 'Position',  [0, 0, 1000, 1000])
% savefig('F_GW_RD');
% saveas(gcf,'F_GW_RD','epsc');
% close all;