function [matrix_RLM, k0] = partial_RD(TFR, n0, q, sigma_s, gamma_TFR)

[Nfft, L] = size(TFR);
matrix_RLM = zeros(size(TFR));

T_size = 150;
Chi2_table = zeros(T_size, 1);
for j=1:T_size
    Chi2_table(j) = chi2inv(1 - 1/(1+2*j), 2);
end

[~, k0] = max(abs(TFR(:, n0)));
rq_n0 = real(q(k0, n0));

    %% local maximum test
    function [LM_bool] = is_LM(k, n)
        if k == Nfft || k == 1
            LM_bool = 0;
            return;
        end
        LM_bool = abs(TFR(k, n)) > abs(TFR(k - 1, n))...
            && abs(TFR(k, n)) > abs(TFR(k + 1, n));
    end

    %% RLM test function
    function [RLM_bool] = is_RLM(k, n, rq)
            std_lcg = 1/(sqrt(2*pi)*sigma_s)*sqrt(1 + sigma_s^4*rq^2);
            std_lcg = ceil(std_lcg*Nfft/L);
            Im = max(1, k -std_lcg):min(Nfft, k +std_lcg);

            TH = gamma_TFR*sqrt(Chi2_table(std_lcg));
            Im_count = sum(abs(TFR(Im, n)) > TH);
            RLM_bool = Im_count > 5;
    end

    %% iteration function
    function partial_RD_iteration(dir, lim)
        % dir = 1 : forward iteration
        % dir = -1 : backward iteration
        LM_k = k0;
        rq = rq_n0;
        for n=(n0 + dir):dir:lim
            %% look for nearest local maximum
            rq_grid = round(rq*Nfft/(L^2));
            next_index = max(1, min(Nfft, LM_k +dir*rq_grid));
            
            V = 1;
            LM_k = next_index;
            while 1
                a = max(1, next_index - V);
                b = min(Nfft, next_index + V);
                if is_LM(b - 1, n)
                    LM_k = b - 1;
                    if is_LM(a + 1, n)
                        if abs(TFR(b - 1)) < abs(TFR(a + 1))
                            LM_k = a + 1;
                        end
                    end
                    break;
                elseif is_LM(a + 1, n)
                    LM_k = a + 1;
                    break;
                elseif a == 1 && b == Nfft
                    % no local maximum in at time index n
                    break;
                end
                V = V + 1;
            end
            
            %% RLM test
            rq = real(q(LM_k, n));
            if is_RLM(LM_k, n, rq)
                matrix_RLM(LM_k, n) = 1;
            else
                break;
            end
        end
    end

if is_RLM(k0, n0, rq_n0)
    matrix_RLM(k0, n0) = 1;
else
    return;
end


%% forward iterations
partial_RD_iteration(1, L);

%% backward iterations
partial_RD_iteration(-1, 1);

end
