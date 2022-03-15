function [global_counts, all_counts] = R1_global_count_rep(s_in, L, Nfft, sigma_s, SNRs, Y_th, N_rep, k_IF_min, k_IF_max)

N_SNRs = length(SNRs);
N_th = length(Y_th);
all_counts = zeros(N_SNRs, N_th, N_rep);
global_counts = zeros(N_SNRs, N_th, N_rep);

for nr = 1:N_rep
    noise = randn(L,1) + 1i*randn(L,1);
    for ns = 1:N_SNRs
        SNR = SNRs(ns);
        s_noise = sigmerge(s_in, noise, SNR);
        [g, Lh] = create_gaussian_window(L, Nfft, sigma_s);
        [STFT] = tfrstft(s_noise, Nfft, 1, g, Lh);
        [STFT_LM] = LM_from_STFT(STFT);
        g_Vg = median(abs(real(STFT(:))))/0.6745;
        
        for ny = 1:N_th
            y = Y_th(ny);
            Cg = y*g_Vg;
            STFT_LM_g = STFT_LM.*(abs(STFT) > Cg);

            count = 0;
            c_global = 0;
            for n=1:L
                k1 = k_IF_min(n);
                k2 = k_IF_max(n);
                [~, k_LM] = max(abs(STFT_LM_g(:, n)));
                if k_LM >= k1 && k_LM <= k2
                    c_global = c_global + 1;
                    count = count + 1;
                    continue;
                end
                c_n = sum(abs(STFT_LM_g(k1:k2, n)));
                if c_n > 0
                    count = count + 1;
                end
            end
            
%             count = SNR;
            all_counts(ns, ny, nr) = count;
            global_counts(ns, ny, nr) = c_global;
        end
    end
    if mod(nr, 10) == 0
        fprintf("%u/%u ", nr, N_rep);
    end
end
fprintf("\n");
end

