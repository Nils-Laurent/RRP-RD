close all;

%% signal definition
L = 4096;
t = (0:L-1)'/L;

B = 3*L/4;
phi_LC = L/8*t+B*(t.^2)/2;
s_LC = exp(2*1i*pi*phi_LC);
sigma_LC = 1/sqrt(B);

phip_LC = L/8 + B*t;

Nfft = 512;
Nr = 1;
cas = 1;
poly_degree = 5;

%% TFR LC
% noise = randn(L,1)+1i*randn(L,1);
% save noise_F1_zero noise
load noise_F1_zero
sn_LC = sigmerge(s_LC, noise, -10);
[g, Lg] = create_gaussian_window(L, Nfft, sigma_LC);

[TFR, ~, ~, q] = FM_operators(sn_LC, L, Nfft, g, Lg, sigma_LC);

S_LM = zeros(size(TFR));

A_TFR = abs(TFR);

LM_TFR = abs(TFR);
for n=1:L
    for k=2:(Nfft-1)
        if A_TFR(k, n) > A_TFR(k - 1, n)...
            && A_TFR(k, n) > A_TFR(k + 1, n)
            % local maximum at (k, n);
            continue;
        end
        
        LM_TFR(k, n) = 0;
    end
    
end
max_LM_v = max(max(LM_TFR));

A = 3;
for n=1:L
    [~, k_vec] = sort(LM_TFR(:, n), 'descend');
    for k=k_vec(A+1:end)
        LM_TFR(k, n) = 0;
    end
end

figure;
imagesc((0:L-1)/L, (0:Nfft-1)*L/Nfft, LM_TFR);
set(gca,'ydir','normal');
colormap(flipud(gray));
axis square;
xlim([0.41 0.89]);
ylim([1800 3250]);
xlabel('time', 'interpreter', 'latex');
ylabel('frequency', 'interpreter', 'latex');
xAX = get(gca,'XAxis');
set(xAX,'FontSize', 26);
yAX = get(gca,'YAxis');
set(yAX,'FontSize', 26);

hold on;
plot((0:L-1)/L, phip_LC, 'r--', 'LineWidth', 2);
hold off;

% set(groot, 'defaultLegendInterpreter', 'latex');
% lgd = legend('Location', 'northwest');
% lgd.FontSize = 24;

pbaspect([1 1 1]);
set(gcf, 'Position',  [0, 0, 1000, 1000])

savefig('F1_R1_zeros_LC');
saveas(gcf,'F1_R1_zeros_LC','epsc');
