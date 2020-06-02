close all;

load fig1-observed-Hanford.txt
t1=fig1_observed_Hanford(:,1);
s1=fig1_observed_Hanford(:,2);

load fig1-waveform-H.txt
% load fig1-residual-H.txt
sw=fig1_waveform_H(:,2)';

t = t1(:);
s = s1(:);
s = hilbert(s);
L = length(s);
sigma_s = 0.05;
T = 0.21;
Nfft = 4096;

[g, Lg] = create_gaussian_window(L, Nfft, sigma_s);

[STFT, ~, ~, QM] = FM_operators(s, Nfft, g, Lg, sigma_s);

IFs_vec = [];
TCs_vec = [];
SNRs = zeros(2, 3);


Nr = 1;
TOLs = [1, 2, 3];
sig_LCR = [];
denoised_STFT = [];
ind = 0;
for TOL = TOLs
    ind = ind + 1;
    TOL
    [Cs, XCs, Qs, TFR_inter] = novel_RRP_RD_spline(STFT(1:128, :), QM(1:128, :), sigma_s, Nr, TOL);
    TCs = (XCs - 1)/L;
    IF_E = (fnval(Qs, (XCs - 1)/L) - 1)*L/Nfft;
    IM_E = fnval(fnder(Qs, 1), (XCs - 1)/L)*L/Nfft; % derivative of Qs
    sigma_vec = 1/(sqrt(2*pi)*sigma_s)*sqrt( 1+sigma_s^4*IM_E.^2 );
    range_vec = ceil(3*sigma_vec*Nfft/L);
    IFs_vec = [IFs_vec; IF_E];
    TCs_vec = [TCs_vec; TCs];

    KY_lower = max(1, Cs - range_vec);
    KY_upper = min(128, Cs + range_vec);
    [m_simple, ~] = MR_simple(STFT(1:128,:), Nfft, XCs, g, Lg, KY_lower, KY_upper, Nr);
    m_simple = real(m_simple);
    [m_LCR, STFT_LCR] = LCR_partial(STFT(1:128,:), g, Lg, Nfft, XCs, IF_E, IM_E, sigma_s);
    m_LCR = real(m_LCR);
    
    if TOL == 3
        sig_LCR = m_LCR;
        denoised_STFT = STFT_LCR;
    end

    
    SNRs(1, ind) = snr(sw(XCs-1), m_simple - sw(XCs-1));
    SNRs(2, ind) = snr(sw(XCs-1), m_LCR - sw(XCs-1));

end

SNRs

%% figures
set(groot, 'defaultLegendInterpreter', 'latex');

figure;
imagesc((0:L-1)/L*T, (0:127)*L/(Nfft*T), abs(STFT(1:128, :)));
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
savefig('F_GW_STFT');
saveas(gcf,'F_GW_STFT','epsc');
close all;

figure;
imagesc((0:L-1)/L*T, (0:127)*L/(Nfft*T), abs(denoised_STFT));
set(gca,'ydir','normal');
colormap(flipud(gray));
axis square
xlabel('time', 'interpreter', 'latex');
ylabel('frequency', 'interpreter', 'latex');
ylim([0, 127*L/(Nfft*T)]);
xAX = get(gca,'XAxis');
set(xAX,'FontSize', 26);
yAX = get(gca,'YAxis');
set(yAX,'FontSize', 26);
pbaspect([1 1 1]);
set(gcf, 'Position',  [0, 0, 1000, 1000])
savefig('F_GW_LCR_STFT');
saveas(gcf,'F_GW_LCR_STFT','epsc');
close all;

figure;
imagesc((0:L-1)/L*T, (0:127)*L/(Nfft*T), abs(TFR_inter));
set(gca,'ydir','normal');
colormap(flipud(gray));
axis square
xlabel('time', 'interpreter', 'latex');
ylabel('frequency', 'interpreter', 'latex');
hold on;
plot(TCs_vec(3, :)*T, IFs_vec(3, :)/T,...
    'LineWidth', 2, 'DisplayName', 'RRP-RD (tol = 3)');

hold off;
lgd = legend;
lgd.FontSize = 37;
xAX = get(gca,'XAxis');
set(xAX,'FontSize', 26);
yAX = get(gca,'YAxis');
set(yAX,'FontSize', 26);
pbaspect([1 1 1]);
set(gcf, 'Position',  [0, 0, 1000, 1000])
savefig('F_GW_RD');
saveas(gcf,'F_GW_RD','epsc');
close all;

figure;
hold on;
plot((1:L-1)/L*T, sw, 'LineWidth', 2,...
    'DisplayName', 'Numerical relativity');
plot(TCs_vec(3, :)*T, sig_LCR, '-.', 'LineWidth', 2,...
    'DisplayName', 'RRP-MR-LCR (tol = 3)');
hold off;
xlabel('time', 'interpreter', 'latex');
ylabel('dB', 'interpreter', 'latex');
xlim([TCs_vec(3, 1)*T, TCs_vec(3, end)*T]);
lgd = legend;
lgd.FontSize = 37;
xAX = get(gca,'XAxis');
set(xAX,'FontSize', 26);
yAX = get(gca,'YAxis');
set(yAX,'FontSize', 26);
pbaspect([1 1 1]);
set(gcf, 'Position',  [0, 0, 1000, 1000])
savefig('F_GW_NR');
saveas(gcf,'F_GW_NR','epsc');
close all;