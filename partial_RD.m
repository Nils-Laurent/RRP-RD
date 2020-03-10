function [Weights, k0] = partial_RD(TFR, n0, q, sigma_s, gamma_TFR, Weights)

[Nfft, L] = size(TFR);

% Values from the Chi-squared table : P[X > Cp] = p
%ratio = 1/100;
%C1 = 9.2103;
% ratio = 1/10;
% C10 = 4.6052;

T_size = 10;
Chi2_table = zeros(T_size, 1);
for j=1:T_size
    Chi2_table(j) = chi2inv(1 - 1/(1+2*j), 2);
end

% TH = gamma_TFR*sqrt(C10);

[~, k0] = max(abs(TFR(:, n0)));
rq_n0 = round(Nfft/(L^2)*real(q(k0, n0)));

    %% local maximum test
    function [LM_bool] = is_LM(k, n)
        LM_bool = abs(TFR(k, n)) > abs(TFR(k - 1, n))...
            && abs(TFR(k, n)) > abs(TFR(k + 1, n));
    end

    %% RLM test function
    function [RLM_bool] = is_RLM(k, n, rq)
            std_lcg = 1/(sqrt(2*pi)*sigma_s)*sqrt(1 + sigma_s^4*rq^2);
            std_lcg = ceil(std_lcg*Nfft/L);
            Im = max(1, k -std_lcg):min(Nfft, k +std_lcg);

            TH = gamma_TFR*sqrt(Chi2_table(std_lcg));
            Im_ratio = sum(abs(TFR(Im, n)) > TH)/length(Im);
            RLM_bool = Im_ratio > 1/(1 + 2*std_lcg);
    end

    %% iteration function
    function partial_RD_iteration(dir, lim)
        % dir = 1 : forward iteration
        % dir = -1 : backward iteration
        k = k0;
        rq = rq_n0;
        for n=(n0 + dir):dir:lim
            %% look for nearest local maximum
            next_index = max(1, min(Nfft, k +dir*rq));
            
            V = 1;
            LM = next_index;
            while 1
                a = max(1, next_index - V);
                b = min(Nfft, next_index + V);
                if is_LM(b - 1, n)
                    LM = b - 1;
                    if is_LM(a + 1, n)
                        if abs(TFR(b - 1)) < abs(TFR(a + 1))
                            LM = a + 1;
                        end
                    end
                    break;
                elseif is_LM(a + 1, n)
                    LM = a + 1;
                    break;
                elseif a == 1 && b == Nfft
                    % no local maximum in at time index n
                    break;
                end
                V = V + 1;
            end
            
            %% RLM test
            rq = round(Nfft/(L^2)*real(q(LM, n)));
            if is_RLM(LM, n, rq)
                Weights(LM, n) = Weights(LM, n) + 1;
            else
                break;
            end
        end
    end

if is_RLM(k0, n0, rq_n0)
    Weights(k0, n0) = Weights(k0, n0) + 1;
else
    return;
end


%% forward iterations
partial_RD_iteration(1, L);

%% backward iterations
partial_RD_iteration(-1, 1);

end
