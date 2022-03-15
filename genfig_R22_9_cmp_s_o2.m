close all;

%% global var
L = 4096;
t = (0:L-1)'/L;

%% Linear Chirp MCS

phi1 = 256*t+2700*(t.^2)/2;
phi2 = 768*t+3000*(t.^2)/2;
IF1 = 256+2700*t;
IF2 = 768+3000*t;
CR1 = 2700*ones(size(t));
CR2 = 3000*ones(size(t));
s1 = exp(2*1i*pi*phi1);
s2 = exp(2*1i*pi*phi2);
modes_LC = [s1, s2].';
IFs_LC = [IF1, IF2].';
CRs_LC = [CR1, CR2].';
sigma_LC = 0.0188;

%% Test
clwin = 10;
Nfft = 512;
% smooth_p = [1 - 10^(-2), 1 - 10^(-3), 1 - 10^(-4)];
smooth_p = [1 - 10^(-4), 1 - 10^(-5), 1 - 10^(-6)];
SNR_IN = -10:-1;
NRep = 40;
% NRep = 1;

modes = modes_LC;

N_SNR = length(SNR_IN);
N_sp = length(smooth_p);
[Nr, L] = size(modes);
s_in = sum(modes, 1);

IFs = IFs_LC;
CRs = CRs_LC;
sigma_s = sigma_LC;

SNR_IF_Cl = zeros(Nr, N_SNR);
SNR_CR_Cl = zeros(Nr, N_SNR);
SNR_IF_RRP = zeros(Nr, N_SNR, 3);
SNR_CR_RRP = zeros(Nr, N_SNR, 3);

% for n=1:N_SNR
%     for k=1:NRep
%         fprintf('snr %d/%d, rep %d/%d\n', n, length(SNR_IN), k, NRep);
%         
%         noise = randn(L, 1)+1i*randn(L, 1);
%         s_noise = sigmerge(transpose(s_in), noise, SNR_IN(n));
%         [g, Lh] = create_gaussian_window(L, Nfft, sigma_s);
%         X_win = 2*Lh:(L-2*Lh);
%         [STFT, omega, omega2, QM, ~, tau] = FM_operators(s_noise, L, Nfft, g, Lh, sigma_s);
% 
%         [Cs_simple] = exridge_mult(STFT, Nr, 0, 0, clwin);
% 
%         %% Classic SNR
%         Cs_IF = zeros(size(Cs_simple));
%         Cs_CR = zeros(size(Cs_simple));
%         for m = 1:Nr
%             ref_mode = modes(m, X_win);
%             ref_IF = IFs(m, X_win);
%             ref_CR = CRs(m, X_win);
%             
%             for n_2=X_win
%                 Cs_IF(m, n_2) = real(omega2(Cs_simple(m, n_2), n_2));
%                 Cs_CR(m, n_2) = real(QM(Cs_simple(m, n_2), n_2));
%                 
%             end
% 
%             x_Cl_IF = snr(ref_IF, Cs_IF(m, X_win) - ref_IF);
%             SNR_IF_Cl(m, n) = SNR_IF_Cl(m, n) + x_Cl_IF/NRep;
%             x_Cl_CR = snr(ref_CR, Cs_CR(m, X_win) - ref_CR);
%             SNR_CR_Cl(m, n) = SNR_CR_Cl(m, n) + x_Cl_CR/NRep;
%         end
% 
%         %% RRP
%         for ns=1:N_sp
%             [Spl, ~] = R1_RRP_RD(STFT, QM, omega, tau, L, Nfft, Nr, sigma_s, smooth_p(ns));
%             
%             for m = 1:Nr
%                 ref_mode = modes(m, X_win);
%                 ref_IF = IFs(m, X_win);
%                 ref_CR = CRs(m, X_win);
%                 
%                 RRP_IF = fnval(Spl(m).spline, (0:L-1)/L);
%                 RRP_CR = fnval(fnder(Spl(m).spline), (0:L-1)/L);
%                 
%                 x_RRP_IF = snr(ref_IF, RRP_IF(X_win) - ref_IF);
%                 SNR_IF_RRP(m, n, ns) = SNR_IF_RRP(m, n, ns) + x_RRP_IF/NRep;
%                 x_RRP_CR = snr(ref_CR, RRP_CR(X_win) - ref_CR);
%                 SNR_CR_RRP(m, n, ns) = SNR_CR_RRP(m, n, ns) + x_RRP_CR/NRep;
%             end
%         end
%     end
% end
% save("data_R22_cmp_s_o2.mat",...
%     'SNR_IF_Cl', 'SNR_CR_Cl', 'SNR_IF_RRP', 'SNR_CR_RRP');
load("data_R22_cmp_s_o2.mat");

%% fig var
B1 = [0 0.4470 0.7410];
B3 = [0.3010 0.7450 0.9330];
B2 = (B3 + B1)/2;

R_MB = [0.6350 0.0780 0.1840];

set(groot, 'defaultLegendInterpreter', 'latex');

m1 = 1;

figure;
hold on;
plot(SNR_IN, SNR_IF_Cl(m1, :),...
    'Linewidth', 2,...
    'DisplayName', '$\hat{\omega}^{[2]}$',...
    'Color', R_MB);
plot(SNR_IN, SNR_IF_RRP(m1, :, 1), '-*',...
    'Linewidth', 2, 'Markersize', 10,...
    'DisplayName', '$s^{fin}$, $\lambda = 10^{-4}$',...
    'Color', B1);
plot(SNR_IN, SNR_IF_RRP(m1, :, 2), '-s',...
    'Linewidth', 2, 'Markersize', 10,...
    'DisplayName', '$s^{fin}$, $\lambda = 10^{-5}$',...
    'Color', B2);
plot(SNR_IN, SNR_IF_RRP(m1, :, 3), '-o',...
    'Linewidth', 2, 'Markersize', 10,...
    'DisplayName', '$s^{fin}$, $\lambda = 10^{-6}$',...
    'Color', B3);
hold off;
lgd = legend('Location', 'southeast');
lgd.FontSize = 36;
xlabel('SNR in', 'interpreter', 'latex');
ylabel('SNR out', 'interpreter', 'latex');
xAX = get(gca,'XAxis');
set(xAX,'FontSize', 26);
yAX = get(gca,'YAxis');
set(yAX,'FontSize', 26);
pbaspect([1 1 1]);
set(gcf, 'Position',  [0, 0, 1000, 1000])
savefig("fig_R22_cmp_s_o2_ifm1");
saveas(gcf, "fig_R22_cmp_s_o2_ifm1", 'epsc');

figure;
hold on;
plot(SNR_IN, SNR_CR_Cl(m1, :),...
    'Linewidth', 2,...
    'DisplayName', '$\hat{q}_{\tilde{f}}$',...
    'Color', R_MB);
plot(SNR_IN, SNR_CR_RRP(m1, :, 1), '-*',...
    'Linewidth', 2, 'Markersize', 10,...
    'DisplayName', "$(s^{fin})\\'$, $\lambda = 10^{-4}$",...
    'Color', B1);
plot(SNR_IN, SNR_CR_RRP(m1, :, 2), '-s',...
    'Linewidth', 2, 'Markersize', 10,...
    'DisplayName', "$(s^{fin})\\'$, $\lambda = 10^{-5}$",...
    'Color', B2);
plot(SNR_IN, SNR_CR_RRP(m1, :, 3), '-o',...
    'Linewidth', 2, 'Markersize', 10,...
    'DisplayName', "$(s^{fin})\\'$, $\lambda = 10^{-6}$",...
    'Color', B3);
hold off;
lgd = legend('Location', 'southeast');
lgd.FontSize = 36;
xlabel('SNR in', 'interpreter', 'latex');
ylabel('SNR out', 'interpreter', 'latex');
xAX = get(gca,'XAxis');
set(xAX,'FontSize', 26);
yAX = get(gca,'YAxis');
set(yAX,'FontSize', 26);
pbaspect([1 1 1]);
set(gcf, 'Position',  [0, 0, 1000, 1000])
savefig("fig_R22_cmp_s_o2_crm1");
saveas(gcf, "fig_R22_cmp_s_o2_crm1", 'epsc');
% return;

m2 = 2;

figure;
hold on;
plot(SNR_IN, SNR_IF_Cl(m2, :),...
    'Linewidth', 2,...
    'DisplayName', '$\hat{\omega}^{[2]}$',...
    'Color', R_MB);
plot(SNR_IN, SNR_IF_RRP(m2, :, 1), '-*',...
    'Linewidth', 2, 'Markersize', 10,...
    'DisplayName', '$s^{fin}$, $\lambda = 10^{-4}$',...
    'Color', B1);
plot(SNR_IN, SNR_IF_RRP(m2, :, 2), '-s',...
    'Linewidth', 2, 'Markersize', 10,...
    'DisplayName', '$s^{fin}$, $\lambda = 10^{-5}$',...
    'Color', B2);
plot(SNR_IN, SNR_IF_RRP(m2, :, 3), '-o',...
    'Linewidth', 2, 'Markersize', 10,...
    'DisplayName', '$s^{fin}$, $\lambda = 10^{-6}$',...
    'Color', B3);
hold off;
lgd = legend('Location', 'southeast');
lgd.FontSize = 36;
xlabel('SNR in', 'interpreter', 'latex');
ylabel('SNR out', 'interpreter', 'latex');
xAX = get(gca,'XAxis');
set(xAX,'FontSize', 26);
yAX = get(gca,'YAxis');
set(yAX,'FontSize', 26);
pbaspect([1 1 1]);
set(gcf, 'Position',  [0, 0, 1000, 1000])
savefig("fig_R22_cmp_s_o2_ifm2");
saveas(gcf, "fig_R22_cmp_s_o2_ifm2", 'epsc');

figure;
hold on;
plot(SNR_IN, SNR_CR_Cl(m2, :),...
    'Linewidth', 2,...
    'DisplayName', '$\hat{q}_{\tilde{f}}$',...
    'Color', R_MB);
plot(SNR_IN, SNR_CR_RRP(m2, :, 1), '-*',...
    'Linewidth', 2, 'Markersize', 10,...
    'DisplayName', "$(s^{fin})\\'$, $\lambda = 10^{-4}$",...
    'Color', B1);
plot(SNR_IN, SNR_CR_RRP(m2, :, 2), '-s',...
    'Linewidth', 2, 'Markersize', 10,...
    'DisplayName', "$(s^{fin})\\'$, $\lambda = 10^{-5}$",...
    'Color', B2);
plot(SNR_IN, SNR_CR_RRP(m2, :, 3), '-o',...
    'Linewidth', 2, 'Markersize', 10,...
    'DisplayName', "$(s^{fin})\\'$, $\lambda = 10^{-6}$",...
    'Color', B3);
hold off;
lgd = legend('Location', 'southeast');
lgd.FontSize = 36;
xlabel('SNR in', 'interpreter', 'latex');
ylabel('SNR out', 'interpreter', 'latex');
xAX = get(gca,'XAxis');
set(xAX,'FontSize', 26);
yAX = get(gca,'YAxis');
set(yAX,'FontSize', 26);
pbaspect([1 1 1]);
set(gcf, 'Position',  [0, 0, 1000, 1000])
savefig("fig_R22_cmp_s_o2_crm2");
saveas(gcf, "fig_R22_cmp_s_o2_crm2", 'epsc');

