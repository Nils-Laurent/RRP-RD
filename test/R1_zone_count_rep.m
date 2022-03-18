function [N_hit_mean, N_hit_var] = R1_zone_count_rep(s_in, L, Nfft, sigma_s, SNRs, N_rep, k_IF_min, k_IF_max)
y = 2;
N_snr = length(SNRs);

N_hit_mean = zeros(3, N_snr);
N_hit_var = zeros(3, N_snr);
for nn = 1:N_snr
    fprintf("snr = %d (%u/%u) : ", SNRs(nn), nn, N_snr);
    tmp_rep = zeros(3, N_rep);
    
    for nr = 1:N_rep
        noise = randn(L,1) + 1i*randn(L,1);
        s_noise = add_noise(s_in, noise, SNRs(nn));

        [STFT, TFR] = sst2(s_noise, sigma_s, Nfft);
        QM = TFR.q_hat;
        omega = TFR.omega1_hat;
        tau = TFR.tau;

        [c_nr] = R1_zone_count(STFT, QM, omega, tau, L, Nfft, y, k_IF_min, k_IF_max);
        
        for p=1:3
            tmp_rep(p, nr) = c_nr(p);
        end
        if mod(nr, 10) == 0
            fprintf("%u/%u ", nr, N_rep);
        end
    end
    fprintf("\n");
    
    for p=1:3
        N_hit_mean(p, nn) = mean(tmp_rep(p, :));
        N_hit_var(p, nn) = std(tmp_rep(p, :));
    end
end

end

