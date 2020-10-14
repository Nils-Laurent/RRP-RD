function [ridge, km] = novel_partial_RD(S_LM, S_LM_sorted, m, km, QM, Nr, sigma_s)

[Nfft, L] = size(S_LM);
ridge = zeros(L, 1);
% matrix_RLM = zeros(size(S_LM));

% [~, km] = max(abs(RLMs(:, m)));
rq_n0 = real(QM(km, m));

    %% iteration function
    function partial_RD_iteration(dir, lim)
        % dir = 1 : forward iteration
        % dir = -1 : backward iteration
        k_next = km;
        rq = rq_n0;
        for n=(m + dir):dir:lim
            %% look for a local maxima close to k_rq
            rq_grid = round(rq*Nfft/(L^2));
            k_rq = max(1, min(Nfft, k_next +dir*rq_grid));
            
            V = 0;
            k_next = k_rq;
            while 1
                a = max(1, k_rq - V);
                b = min(Nfft, k_rq + V);
                if S_LM(b, n) > 0
                    k_next = b;
                    if S_LM(a, n) > 0
                        if S_LM(b, n) < S_LM(a, n)
                            k_next = a;
                        end
                    end
                    break;
                elseif S_LM(a, n) > 0
                    k_next = a;
                    break;
                elseif a == 1 && b == Nfft
                    % no local maximum in at time index n
                    k_next = 0;
                    break;
                end
                V = V + 1;
            end
            
            %% test k_next
            if k_next == 0
                break;
            end
            
            A = 2*Nr + 1;
            %[~, k_vec] = sort(S_LM(:, n), 'descend');
            k_vec = S_LM_sorted(:, n);
            k_ = 0;
            for p=1:A
                if k_vec(p) == k_next
                    k_ = 1;
                end
            end
            
            if k_ == 0
                break;
            end
            
%             std_lc_g = 1/(sqrt(2*pi)*sigma_s)*sqrt(1 + sigma_s^4*rq^2);
%             std_lc_g = ceil(std_lc_g*Nfft/L);
%             if abs(next_index - RLM_k) > std_lc_g
%                 break;
%             end
            
            
            rq = real(QM(k_next, n));
            ridge(n) = k_next;
        end
    end

if S_LM(km, m) > 0
%     matrix_RLM(km, m) = 1;
    ridge(m) = km;
else
    return;
end


%% forward iterations
partial_RD_iteration(1, L);

%% backward iterations
partial_RD_iteration(-1, 1);

end
