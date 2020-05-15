close all;

%% signal definition
L = 4096;
t = (0:L-1)'/L;

B = 3*L/4;
phi_LC = L/8*t+B*(t.^2)/2;
s_LC = exp(2*1i*pi*phi_LC);
sigma_LC = 1/sqrt(B);

B = floor(2*L/(2*pi));
phi_cos = L/2*t+B/(2*pi)*cos(2*pi*t);
s_cos = exp(2*1i*pi*phi_cos);
sigma_cos = 0.0142;

B = log(4096 - 510);
phi_exp = 500*t+exp(B*t)/B;
s_exp = exp(2*1i*pi*phi_exp);
sigma_exp = 0.028;

Nfft = 512;
Nr = 1;
cas = 1;
poly_degree = 5;

%% TFR LC
% noise = randn(L,1)+1i*randn(L,1);
% save noise_SC_LC noise
load noise_SC_cos
sn_LC = sigmerge(s_LC, noise, -10);
[g, Lg] = create_gaussian_window(L, Nfft, sigma_LC);

[TFR, ~, ~, q] = FM_operators(sn_LC, Nfft, g, Lg, sigma_LC);
[Cs, Qs, KY_lower, KY_upper, ~] = novel_RRP_RD(TFR, q, sigma_LC, Nr, poly_degree);

fig_TFR(TFR);
savefig('F2_TFR_LC');
saveas(gcf,'F2_TFR_LC','epsc');
close all;

fig_RD(TFR, Cs, Qs, KY_lower, KY_upper);
savefig('F2_RD_LC');
saveas(gcf,'F2_RD_LC','epsc');
close all;

%% TFR cos
% noise = randn(L,1)+1i*randn(L,1);
% save noise_SC_cos noise
load noise_SC_cos
sn_cos = sigmerge(s_cos, noise, -10);
[g, Lg] = create_gaussian_window(L, Nfft, sigma_cos);

[TFR, ~, ~, q] = FM_operators(sn_cos, Nfft, g, Lg, sigma_cos);
[Cs, Qs, KY_lower, KY_upper, ~] = novel_RRP_RD(TFR, q, sigma_cos, Nr, poly_degree);

fig_TFR(TFR);
savefig('F2_TFR_cos');
saveas(gcf,'F2_TFR_cos','epsc');
close all;

fig_RD(TFR, Cs, Qs, KY_lower, KY_upper);
savefig('F2_RD_cos');
saveas(gcf,'F2_RD_cos','epsc');
close all;

%% TFR exp
% noise = randn(L,1)+1i*randn(L,1);
% save noise_SC_exp noise
load noise_SC_exp
sn_exp = sigmerge(s_exp, noise, -10);
[g, Lg] = create_gaussian_window(L, Nfft, sigma_exp);

[TFR, ~, ~, q] = FM_operators(sn_exp, Nfft, g, Lg, sigma_exp);
[Cs, Qs, KY_lower, KY_upper, ~] = novel_RRP_RD(TFR, q, sigma_exp, Nr, poly_degree);

fig_TFR(TFR);
savefig('F2_TFR_exp');
saveas(gcf,'F2_TFR_exp','epsc');
close all;

fig_RD(TFR, Cs, Qs, KY_lower, KY_upper);
savefig('F2_RD_exp');
saveas(gcf,'F2_RD_exp','epsc');
close all;

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
    
    Cs(Cs == 0) = 1;
    Cs = Cs - 1;
    Cs = Cs*L/Nfft;
    Cs(Cs == 0) = nan;
    
    QY = polyval(Qs(1, :), (0:L-1)/L);
    QY_lower = KY_lower*L/Nfft;
    QY_upper = KY_upper*L/Nfft;
    
    figure;
    imagesc((0:L-1)/L, (0:Nfft-1)*L/Nfft, abs(TFR));
    set(gca,'ydir','normal');
    colormap(flipud(gray));
    hold on;
    plot((0:L-1)/L, QY_lower, 'r-.',...
        'DisplayName', '$TF^-$', 'LineWidth', 1);
    plot((0:L-1)/L, QY_upper, 'r--',...
        'DisplayName', '$TF^+$', 'LineWidth', 1);
    plot((0:L-1)/L, QY, 'r',...
        'DisplayName', '$D^{fin}$', 'LineWidth', 1);
    plot((0:L-1)/L, Cs, 'c',...
        'DisplayName', '$\mathcal{M}^*$', 'LineWidth', 2);
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
