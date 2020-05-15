close all;

%% signal definition
L = 4096;
t = (0:L-1)'/L;

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

B1 = log(2000);
B2 = log(2800);
phi1 = 2000*t+exp(B1*t)/B1;
phi2 = 200*t+exp(B2*t)/B2;
s1 = exp(2*1i*pi*phi1);
s2 = exp(2*1i*pi*phi2);
s_exp = s1 + s2;
sigma_exp = 0.0336;

Nr = 2;
Nfft = 512;
cas = 1;
poly_degree = 5;

SNR_IN = -10;

%% TFR LC

% FILE_ = load("noise_poly");
% noise = FILE_.noise;

% noise = randn(L,1)+1i*randn(L,1);
% save noise_MCS_LC noise
load noise_MCS_LC
s_noise = sigmerge(s_LC, noise, SNR_IN);
[g, Lg] = create_gaussian_window(L, Nfft, sigma_LC);

[TFR, ~, ~, q] = FM_operators(s_noise, Nfft, g, Lg, sigma_LC);
[Cs, Qs, KY_lower, KY_upper, ~] = novel_RRP_RD(TFR, q, sigma_LC, Nr, poly_degree);

fig_TFR(TFR);
savefig('F3_TFR_LC');
saveas(gcf,'F3_TFR_LC','epsc');
close all;

fig_RD(TFR, Cs, Qs, KY_lower, KY_upper);
savefig('F3_RD_LC');
saveas(gcf,'F3_RD_LC','epsc');
close all;

%% TFR cos lim
% noise = randn(L,1)+1i*randn(L,1);
% save noise_MCS_cos_lim noise
load noise_MCS_cos_lim
s_noise = sigmerge(s_cos_lim, noise, SNR_IN);
[g, Lg] = create_gaussian_window(L, Nfft, sigma_cos_lim);

[TFR, ~, ~, q] = FM_operators(s_noise, Nfft, g, Lg, sigma_cos_lim);
[Cs, Qs, KY_lower, KY_upper, ~] = novel_RRP_RD(TFR, q, sigma_cos_lim, Nr, poly_degree);

fig_TFR(TFR);
savefig('F3_TFR_cos_lim');
saveas(gcf,'F3_TFR_cos_lim','epsc');
close all;

fig_RD(TFR, Cs, Qs, KY_lower, KY_upper);
savefig('F3_RD_cos_lim');
saveas(gcf,'F3_RD_cos_lim','epsc');
close all;

%% TFR exp LC
% noise = randn(L,1)+1i*randn(L,1);
% save noise_MCS_exp_LC noise
load noise_MCS_exp_LC
s_noise = sigmerge(s_exp_LC, noise, SNR_IN);
[g, Lg] = create_gaussian_window(L, Nfft, sigma_exp_LC);

[TFR, ~, ~, q] = FM_operators(s_noise, Nfft, g, Lg, sigma_exp_LC);
[Cs, Qs, KY_lower, KY_upper, ~] = novel_RRP_RD(TFR, q, sigma_exp_LC, Nr, poly_degree);

fig_TFR(TFR);
savefig('F3_TFR_exp_LC');
saveas(gcf,'F3_TFR_exp_LC','epsc');
close all;

fig_RD(TFR, Cs, Qs, KY_lower, KY_upper);
savefig('F3_RD_exp_LC');
saveas(gcf,'F3_RD_exp_LC','epsc');
close all;

%% TFR exp
% noise = randn(L,1)+1i*randn(L,1);
% save noise_MCS_exp noise
load noise_MCS_exp
s_noise = sigmerge(s_exp, noise, SNR_IN);
[g, Lg] = create_gaussian_window(L, Nfft, sigma_exp);

[TFR, ~, ~, q] = FM_operators(s_noise, Nfft, g, Lg, sigma_exp);
[Cs, Qs, KY_lower, KY_upper, ~] = novel_RRP_RD(TFR, q, sigma_exp, Nr, poly_degree);

fig_TFR(TFR);
savefig('F3_TFR_exp');
saveas(gcf,'F3_TFR_exp','epsc');
close all;

fig_RD(TFR, Cs, Qs, KY_lower, KY_upper);
savefig('F3_RD_exp');
saveas(gcf,'F3_RD_exp','epsc');
% close all;

%% fig functions
function fig_TFR(TFR)
    set(groot, 'defaultLegendInterpreter', 'latex');
    [Nfft, L] = size(TFR);
    
    figure;
    imagesc((0:L-1)/L, (0:Nfft-1)*L/Nfft, abs(TFR));
    set(gca,'ydir','normal');
    colormap(flipud(gray));
    axis square
    xlabel('time', 'interpreter', 'latex');
    ylabel('frequency', 'interpreter', 'latex');
    xAX = get(gca,'XAxis');
    set(xAX,'FontSize', 26);
    yAX = get(gca,'YAxis');
    set(yAX,'FontSize', 26);
    pbaspect([1 1 1]);
    set(gcf, 'Position',  [0, 0, 1000, 1000])
end

function fig_RD(TFR, Cs, Qs, KY_lower, KY_upper)
    set(groot, 'defaultLegendInterpreter', 'latex');
    [Nfft, L] = size(TFR);
    time_vec = (0:(L-1))/L;
    
    Cs(Cs == 0) = 1;
    Cs = Cs - 1;
    Cs = Cs*L/Nfft;
    Cs(Cs == 0) = nan;
    
    QY1 = polyval(Qs(1, :), time_vec);
    QY2 = polyval(Qs(2, :), time_vec);
    QY_lower = KY_lower*L/Nfft;
    QY_upper = KY_upper*L/Nfft;
    
    figure;
    imagesc(time_vec, (0:(Nfft-1))*L/Nfft, abs(TFR));
    set(gca,'ydir','normal');
    colormap(flipud(gray));
    
    hold on;
    plot(time_vec, QY_lower(1, :), 'r-.',...
        'DisplayName', '$TF^-_1$', 'LineWidth', 1);
    plot(time_vec, QY_upper(1, :), 'r--',...
        'DisplayName', '$TF^+_2$', 'LineWidth', 1);
    plot(time_vec, QY1, 'r',...
        'DisplayName', '$D_1^{fin}$', 'LineWidth', 1);
    plot(time_vec, Cs(1, :), 'c',...
        'DisplayName', '$\mathcal{M}^*_1$', 'LineWidth', 2);
    
    plot(time_vec, QY_lower(2, :), 'b-.',...
        'DisplayName', '$TF^-_2$', 'LineWidth', 1);
    plot(time_vec, QY_upper(2, :), 'b--',...
        'DisplayName', '$TF^+_2$', 'LineWidth', 1);
    plot(time_vec, QY2, 'b',...
        'DisplayName', '$D_2^{fin}$', 'LineWidth', 1);
    plot(time_vec, Cs(2, :), 'y',...
        'DisplayName', '$\mathcal{M}^*_2$', 'LineWidth', 2);
    hold off;
    
    xlabel('time', 'interpreter', 'latex');
    ylabel('frequency', 'interpreter', 'latex');
    lgd = legend('Location', 'northwest');
    lgd.FontSize = 24;
    xAX = get(gca,'XAxis');
    set(xAX,'FontSize', 26);
    yAX = get(gca,'YAxis');
    set(yAX,'FontSize', 26);
    pbaspect([1 1 1]);
    set(gcf, 'Position',  [0, 0, 1000, 1000])
end