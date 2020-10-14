function [S_LM] = novel_LM_sum(TFR, QM, sigma_s)

[Nfft, L] = size(TFR);

% fprintf("compute local maxima energies\n");
% gamma_TFR = median(abs(real(TFR(:))))/0.6745;

% T_size = 150;
% Chi2_table = zeros(T_size, 1);
% for j=1:T_size
%     Chi2_table(j) = chi2inv(1 - 1/(1+2*j), 2);
% end

function [LM_bool] = is_LM(k, n)
    if k == Nfft || k == 1
        LM_bool = 0;
        return;
    end
    LM_bool = abs(TFR(k, n)) > abs(TFR(k - 1, n))...
        && abs(TFR(k, n)) > abs(TFR(k + 1, n));
end
function [Im_sum] = get_r(k, n, rq)
        std_lc_g = 1/(sqrt(2*pi)*sigma_s)*sqrt(1 + sigma_s^4*rq^2);
        std_lc_g = ceil(std_lc_g*Nfft/L);
        Im = max(1, k -std_lc_g):min(Nfft, k +std_lc_g);

%         if std_lc_g <= T_size
%             TH = gamma_TFR*sqrt(Chi2_table(std_lc_g));
%         else
%             TH = gamma_TFR*sqrt(chi2inv(1 - 1/(1+2*std_lc_g), 2));
%         end
%         Im_sum = sum((abs(TFR(Im, n)).^2).*(abs(TFR(Im, n)) > TH));
        Im_sum = sum(abs(TFR(Im, n)).^2);
end

S_LM = zeros(size(TFR));

for n=1:L
    for k=2:(Nfft-1)
        if is_LM(k, n)
            S_LM(k, n) = get_r(k, n, real(QM(k, n)));
        end
    end
end

end

