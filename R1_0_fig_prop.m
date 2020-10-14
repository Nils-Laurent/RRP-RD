close all;

%% signal definition
L = 4096;
t = (0:L-1)'/L;

B = floor(2*L/(2*pi));
phi_s = L/2*t+B/(2*pi)*cos(2*pi*t);
s_in = exp(2*1i*pi*phi_s);
sigma_s = 0.0142;

phip_s = L/2-B*sin(2*pi*t);
phipp_s = -2*pi*cos(2*pi*t);

std_s = 1/(sqrt(2*pi)*sigma_s)*sqrt(1 + sigma_s^4*phipp_s.^2);

%% STFT
Nfft = 512;
Nr = 1;
cas = 1;

fx = (0:Nfft-1)*L/Nfft;

k_IF_min = round((phip_s - std_s)*Nfft/L) + 1;
k_IF_max = round((phip_s + std_s)*Nfft/L) + 1;

%% Loops

SNRs = -10:1:2;
N_snr = length(SNRs);
N_rep = 30;
% N_rep = 1;

y = 2.35;
% y0 = 2.35;
% y1 = 3;
% 
% utils = [2536.06666666667 2975.46666666667 3357.26666666667 3634.10000000000...
%     3829.80000000000 3942.16666666667 4001.33333333333 4035 4054.53333333333...
%     4066.20000000000 4068.96666666667];
% utils = (utils - utils(1))/(utils(end) - utils(1));
% Y_vec = y0 + utils*(y1 - y0);
% N_Y = length(Y_vec);

% figure;
% plot(Y_vec);
% pause;
% return;

N_hit_mean = zeros(3, N_snr);
N_hit_var = zeros(3, N_snr);

%% test max
% R1_TF_in_STD(s_in, L, Nfft, sigma_s, y, -10, k_IF_min, k_IF_max)

for nn = 1:N_snr
    fprintf("snr = %d (%u/%u) : ", SNRs(nn), nn, N_snr);
    tmp_rep = zeros(3, N_rep);
    
    for nr = 1:N_rep
        [c_nr] = R1_TF_in_STD(s_in, L, Nfft, sigma_s, y, SNRs(nn), k_IF_min, k_IF_max);
        for p=1:3
            tmp_rep(p, nr) = c_nr(p);
        end
        if mod(nr, 10) == 0
            fprintf("%u/%u ", nr, N_rep);
        end
    end
    fprintf("\n");
    
    for p=1:3
        N_hit_mean(p, nn) = mean(tmp_rep(p, :));
        N_hit_var(p, nn) = std(tmp_rep(p, :));
    end
end

%% figures

figure;
hold on;
plot(SNRs, N_hit_mean(1, :), 'LineWidth', 2,...
    'DisplayName', 'Mean-Max');
plot(SNRs, N_hit_mean(2, :), 'y.-', 'LineWidth', 2,...
    'DisplayName', 'Mean-RP');
plot(SNRs, N_hit_mean(3, :), 'r--', 'LineWidth', 2,...
    'DisplayName', 'Mean-Basins');
hold off;
xlabel('SNRs', 'interpreter', 'latex');
ylabel('mean', 'interpreter', 'latex');
xAX = get(gca,'XAxis');
set(xAX,'FontSize', 26);
yAX = get(gca,'YAxis');
set(yAX,'FontSize', 26);
set(groot, 'defaultLegendInterpreter', 'latex');
lgd = legend('Location', 'southeast');
lgd.FontSize = 24;
pbaspect([1 1 1]);
set(gcf, 'Position',  [0, 0, 1000, 1000]);
% savefig('fig_R1_prop_mean');
% saveas(gcf,'F1_zeros_LC','epsc');
% close all

figure;
hold on;
plot(SNRs, N_hit_var(1, :), 'LineWidth', 2,...
    'DisplayName', 'Std-Max');
plot(SNRs, N_hit_var(2, :), 'y.-', 'LineWidth', 2,...
    'DisplayName', 'Std-RP');
plot(SNRs, N_hit_var(3, :), 'r--', 'LineWidth', 2,...
    'DisplayName', 'Std-Basins');
hold off;
xlabel('SNRs', 'interpreter', 'latex');
ylabel('std', 'interpreter', 'latex');
xAX = get(gca,'XAxis');
set(xAX,'FontSize', 26);
yAX = get(gca,'YAxis');
set(yAX,'FontSize', 26);
set(groot, 'defaultLegendInterpreter', 'latex');
lgd = legend('Location', 'southeast');
lgd.FontSize = 24;
pbaspect([1 1 1]);
set(gcf, 'Position',  [0, 0, 1000, 1000]);
% savefig('fig_R1_prop_var');
% saveas(gcf,'F1_zeros_LC','epsc');
% close all
