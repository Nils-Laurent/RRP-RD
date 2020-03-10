function [Cs] = exridge_n2_MCS(TFR, q, Nr, sigma_s)

[Nfft, L] = size(TFR);
gamma_TFR = median(abs(real(TFR(:))))/0.6745;

Cs_tmp = zeros(Nr, L);

TFR_r = TFR;
for r=1:Nr
    [C_opt_r] = exridge_n2(TFR_r, q, sigma_s, gamma_TFR);
    Cs_tmp(r, :) = C_opt_r;
    
    for n=1:L
        Ck = C_opt_r(n);
        if Ck == 0
            continue;
        end
        rq = round(Nfft/(L^2)*real(q(Ck, n)));
        eta_lim = round(3/(sqrt(2*pi)*sigma_s)*sqrt(1 + sigma_s^4*rq^2)*Nfft/L);
        R_top = min(Nfft, Ck + eta_lim);
        R_bot = max(1, Ck - eta_lim);
        TFR_r(R_bot:R_top, n) = 0;
    end
    
%     XX = Cs_tmp(1, :);
%     XX(XX == 0) = nan;
%     figure;
%     imagesc(1:L, 1:Nfft, abs(TFR_r));
%     set(gca,'ydir','normal');
%     axis square
%     colormap(flipud(gray));
%     hold on;
%     plot(1:L, XX);
%     hold off;
%     pause;
end

% XX = Cs_tmp;
% XX(XX == 0) = nan;
% figure;
% imagesc(1:L, 1:Nfft, abs(TFR));
% set(gca,'ydir','normal');
% axis square
% colormap(flipud(gray));
% hold on;
% plot(1:L, XX);
% hold off;
% pause;

Cs = zeros(Nr, L);
for n=1:L
    Cs(:, n) = sort(Cs_tmp(:, n));
end

end