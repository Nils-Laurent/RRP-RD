close all;

%% signal definition
L = 4096;
t = (0:L-1)'/L;

A = L/32;
B = L - 2*A;

phi_s = A*t+B/2*(t.^2);
s_in = exp(2*1i*pi*phi_s);
sigma_s = 1/sqrt(B);

phip_s = A + B*t;
phipp_s = B*ones(L, 1);

std_s = 1/(sqrt(2*pi)*sigma_s)*sqrt(1 + sigma_s^4*phipp_s.^2);

%% STFT
Nfft = 512;
Nr = 1;
cas = 1;

fx = (0:Nfft-1)*L/Nfft;

k_IF_min = round((phip_s - std_s)*Nfft/L) + 1;
k_IF_max = round((phip_s + std_s)*Nfft/L) + 1;

% [g, Lh] = create_gaussian_window(L, Nfft, sigma_s);
% [STFT] = tfrstft(s_in, Nfft, 1, g, Lh);
% figure;
% imagesc(abs(STFT));
% set(gca,'ydir','normal');
% colormap(flipud(gray));
% axis square;
% hold on;
% plot(k_IF_min, 'r');
% plot(round(phip_s*Nfft/L) + 1);
% plot(k_IF_max, 'r');
% hold off;
% return;

%% Loops

SNR = -10;
Y_th = 2:0.1:3;
N_th = length(Y_th);
N_rep = 100;

N_hit_mean = zeros(3, N_th);
N_hit_var = zeros(3, N_th);

%%
% R1_TF_in_STD(s_in, L, Nfft, sigma_s, 3, SNR, k_IF_min, k_IF_max)

for ny = 1:N_th
    fprintf("y = %d (%u/%u) : ", Y_th(ny), ny, N_th);
    tmp_rep = zeros(3, N_rep);
    for nr = 1:N_rep
        [c_nr] = R1_TF_in_STD(s_in, L, Nfft, sigma_s, Y_th(ny), SNR, k_IF_min, k_IF_max);
        for p=1:3
            tmp_rep(p, nr) = c_nr(p);
        end
        if mod(nr, 10) == 0
            fprintf("%u/%u ", nr, N_rep);
        end
%         pause;
    end
    fprintf("\n");
    
    for p=1:3
        N_hit_mean(p, ny) = mean(tmp_rep(p, :));
        N_hit_var(p, ny) = std(tmp_rep(p, :));
    end
end

figure;
hold on;
plot(Y_th, N_hit_mean(1, :), 'LineWidth', 2,...
    'DisplayName', 'Mean-Max');
plot(Y_th, N_hit_mean(2, :), 'y.-', 'LineWidth', 2,...
    'DisplayName', 'Mean-RP');
plot(Y_th, N_hit_mean(3, :), 'r--', 'LineWidth', 2,...
    'DisplayName', 'Mean-Basins');
hold off;
xlabel('y', 'interpreter', 'latex');
ylabel('mean', 'interpreter', 'latex');
xAX = get(gca,'XAxis');
set(xAX,'FontSize', 26);
yAX = get(gca,'YAxis');
set(yAX,'FontSize', 26);
set(groot, 'defaultLegendInterpreter', 'latex');
lgd = legend('Location', 'northeast');
lgd.FontSize = 24;
pbaspect([1 1 1]);
set(gcf, 'Position',  [0, 0, 1000, 1000])

figure;
hold on;
plot(Y_th, N_hit_var(1, :), 'LineWidth', 2,...
    'DisplayName', 'Std-Max');
plot(Y_th, N_hit_var(2, :), 'y.-', 'LineWidth', 2,...
    'DisplayName', 'Std-RP');
plot(Y_th, N_hit_var(3, :), 'r--', 'LineWidth', 2,...
    'DisplayName', 'Std-Basins');
hold off;
xlabel('y', 'interpreter', 'latex');
ylabel('std', 'interpreter', 'latex');
xAX = get(gca,'XAxis');
set(xAX,'FontSize', 26);
yAX = get(gca,'YAxis');
set(yAX,'FontSize', 26);
set(groot, 'defaultLegendInterpreter', 'latex');
lgd = legend('Location', 'northeast');
lgd.FontSize = 24;
pbaspect([1 1 1]);
set(gcf, 'Position',  [0, 0, 1000, 1000])

