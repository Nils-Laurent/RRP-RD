close all;

%% global var

L = 4096;
t = (0:L-1)'/L;

%% Signal def

B1 = log(2000);
phi1 = 200*t+2300*(t.^2)/2;
phi2 = 2000*t+exp(B1*t)/B1;
s1 = exp(2*1i*pi*phi1);
s2 = exp(2*1i*pi*phi2);
s_exp = s1 + s2;
sigma_exp = 0.0241;

%% test
SNR_in = -10;
Nfft = 512;
smooth_p = 1 - 10^(-4);
Nr = 2;

%% exp
load('noise_R1_TFR_RD_exp.mat');
s_noise = sigmerge(s_exp, noise, SNR_in);
[g, Lh] = create_gaussian_window(L, Nfft, sigma_exp);
[STFT, omega, ~, QM, ~, tau] = FM_operators(s_noise, L, Nfft, g, Lh, sigma_exp);
[Spl_exp, ~] = R1_RRP_RD(STFT, QM, omega, tau, L, Nfft, Nr, sigma_exp, smooth_p);

Spl = Spl_exp;
sigma_s = sigma_exp;
IF1 = fnval(Spl(1).spline, t);
DF1 = fnval(fnder(Spl(1).spline), t);
R1 = 1/(sqrt(2*pi)*sigma_s)*sqrt(1 + sigma_s^4*DF1.^2);
IF2 = fnval(Spl(2).spline, t);
DF2 = fnval(fnder(Spl(2).spline), t);
R2 = 1/(sqrt(2*pi)*sigma_s)*sqrt(1 + sigma_s^4*DF2.^2);


set(groot, 'defaultLegendInterpreter', 'latex');

figure;
imagesc(t, (0:Nfft-1)*L/Nfft, abs(STFT));
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
plot(t, IF1, 'k',...
    'Linewidth', 3,...
    'DisplayName', '$s_1^{fin}$');
plot(t, IF2, 'k--',...
    'Linewidth', 3,...
    'DisplayName', '$s_2^{fin}$');
hold off;
lgd = legend('Location', 'southeast');
lgd.FontSize = 24;
pbaspect([1 1 1]);
set(gcf, 'Position',  [0, 0, 1000, 1000])