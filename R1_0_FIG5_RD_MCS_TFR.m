close all;

%% global var

L = 4096;
t = (0:L-1)'/L;

%% Signal def
% B = 3*L/4;
% phi_LC = L/8*t+B*(t.^2)/2;
% s_LC = exp(2*1i*pi*phi_LC);
% sigma_LC = 1/sqrt(B);
% 
% B = floor(2*L/(2*pi));
% phi_cos = L/2*t+B/(2*pi)*cos(2*pi*t);
% s_cos = exp(2*1i*pi*phi_cos);
% sigma_cos = 0.0142;
% 
% B = log(4096 - 510);
% phi_exp = 500*t+exp(B*t)/B;
% s_exp = exp(2*1i*pi*phi_exp);
% sigma_exp = 0.028;

phi1 = 256*t+2700*(t.^2)/2;
phi2 = 768*t+3000*(t.^2)/2;
s1 = exp(2*1i*pi*phi1);
s2 = exp(2*1i*pi*phi2);
s_LC = s1 + s2;
sigma_LC = 0.0188;

phi1 = 1400*t+1350/(2*pi)*cos(2*pi*t + pi/2);
phi2 = 3400*t+550/(2*pi)*cos(2*pi*t);
s1 = exp(2*1i*pi*phi1);
s2 = exp(2*1i*pi*phi2);
s_cos = s1 + s2;
sigma_cos = 0.0175;

B1 = log(2000);
phi1 = 200*t+2300*(t.^2)/2;
phi2 = 2000*t+exp(B1*t)/B1;
s1 = exp(2*1i*pi*phi1);
s2 = exp(2*1i*pi*phi2);
s_exp = s1 + s2;
sigma_exp = 0.0241;

%% test
SNR_in = -8;
Nfft = 512;
smooth_p = 1 - 10^(-2);
Nr = 2;

% 3 FIG sig. MCS
% 3 Courbes (spline, std sup, std inf)

%% LC

%% cos

%% exp
noise = randn(L, 1)+1i*randn(L, 1);
s_noise = sigmerge(s_exp, noise, SNR_in);
[g, Lh] = create_gaussian_window(L, Nfft, sigma_exp);
[STFT, omega, ~, QM, ~, tau] = FM_operators(s_noise, L, Nfft, g, Lh, sigma_exp);
[Spl_exp, ~] = R1_RRP_RD(STFT, QM, omega, tau, L, Nfft, Nr, sigma_exp, smooth_p);

IF1 = fnval(Spl_exp(1).spline, t);
DF1 = fnval(fnder(Spl_exp(1).spline), t);
R1 = 1/(sqrt(2*pi)*sigma_LC)*sqrt(1 + sigma_LC^4*DF1.^2);
IF2 = fnval(Spl_exp(2).spline, t);
DF2 = fnval(fnder(Spl_exp(2).spline), t);
R2 = 1/(sqrt(2*pi)*sigma_LC)*sqrt(1 + sigma_LC^4*DF2.^2);

figure;
imagesc((0:L-1)/L, (0:Nfft-1)*L/Nfft, abs(STFT));
set(gca,'ydir','normal');
colormap(flipud(gray));
axis square;
xlabel('time', 'interpreter', 'latex');
ylabel('frequency', 'interpreter', 'latex');
xAX = get(gca,'XAxis');
set(xAX,'FontSize', 26);
yAX = get(gca,'YAxis');
set(yAX,'FontSize', 26);
hold on;
Y = [0.9290 0.6940 0.1250];
plot(t, IF1 - R1, 'b-',...
    'Linewidth', 2,...
    'DisplayName', '$F_1^-$');
plot(t, IF1, 'c',...
    'Linewidth', 2,...
    'DisplayName', '$D_1^{fin}$');
plot(t, IF1 + R1, 'r-',...
    'Linewidth', 2,...
    'DisplayName', '$F_1^+$');
plot(t, IF2 - R2, 'b--',...
    'Linewidth', 2,...
    'DisplayName', '$F_2^-$');
plot(t, IF2, 'c--',...
    'Linewidth', 2,...
    'DisplayName', '$D_2^{fin}$');
plot(t, IF2 + R2, 'r--',...
    'Linewidth', 2,...
    'DisplayName', '$F_2^+$');
hold off;
lgd = legend('Location', 'southeast');
lgd.FontSize = 24;
xlabel('time', 'interpreter', 'latex');
ylabel('frequency', 'interpreter', 'latex');
xAX = get(gca,'XAxis');
set(xAX,'FontSize', 26);
yAX = get(gca,'YAxis');
set(yAX,'FontSize', 26);
pbaspect([1 1 1]);
set(gcf, 'Position',  [0, 0, 1000, 1000])