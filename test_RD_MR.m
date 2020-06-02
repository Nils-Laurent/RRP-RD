function [SNR_NEW_LCR, SNR_NEW_s, SNR_S_RD, SNR_MB_RD,...
    SNR_IF_NEW, SNR_IF_S_RD, SNR_IF_MB_RD] = test_RD_MR(...
    modes, IFs, clwin, sigma_s, Nfft, poly_d, SNR_IN, NRep)

cas = 1;
[Nr, L] = size(modes);
s_in = sum(modes, 1);

SNR_NEW_LCR = zeros(Nr, length(SNR_IN));
SNR_NEW_s = zeros(Nr, length(SNR_IN));
SNR_S_RD = zeros(Nr, length(SNR_IN));
SNR_MB_RD = zeros(Nr, length(SNR_IN));

SNR_IF_NEW = zeros(Nr, length(SNR_IN));
SNR_IF_S_RD = zeros(Nr, length(SNR_IN));
SNR_IF_MB_RD = zeros(Nr, length(SNR_IN));

for k=1:length(SNR_IN)
    for l=1:NRep
        fprintf('snr %d/%d, rep %d/%d\n', k, length(SNR_IN), l, NRep);
        
        noise = randn(L, 1)+1i*randn(L, 1);
        s_noise = sigmerge(transpose(s_in), noise, SNR_IN(k));
        [g, Lg] = create_gaussian_window(L, Nfft, sigma_s);
        [STFT, ~, ~, QM] = FM_operators(s_noise, Nfft, g, Lg, sigma_s);

        fprintf('VFB MB, ');
        [Cs_VFB_MB] = VFB_MB_exridge_MCS(STFT, sigma_s, QM, 2, Nr);
        [K_lower, K_upper] = MR_interval(Cs_VFB_MB, QM, Nfft, sigma_s);
        [m_MB_RD, ~] = MR_simple(STFT, Nfft, 1:L, g, Lg, K_lower, K_upper, Nr);

        fprintf('classic, ');
        [Cs_simple] = exridge_mult(STFT, Nr, 0, 0, clwin);
        [K_lower, K_upper] = MR_interval(Cs_simple, QM, Nfft, sigma_s);
        [m_C_RD, ~] = MR_simple(STFT, Nfft, 1:L, g, Lg, K_lower, K_upper, Nr);

        fprintf('RRP RD, ');
        [~, Qs, KY_lower, KY_upper, ~] = novel_RRP_RD(STFT, QM, sigma_s, Nr, poly_d);
        Cs = zeros(size(modes));
        IF_E = zeros(size(modes));
        IM_E = zeros(size(modes));
        for p = 1:Nr
            PY = polyval(Qs(p, :), (0:L-1)/L);
            IF_E(p, :) = PY;
            IM_E(p, :) = polyval(polyder(Qs(p, :)), (0:L-1)/L);
            Cs(p, :) = round(PY*Nfft/L) + 1;
            Cs(p, :) = min(Nfft, max(1, Cs(p, :)));
        end

        fprintf("LCR\n");
        [m_NEW_simple, ~] = MR_simple(STFT, Nfft, 1:L, g, Lg, KY_lower, KY_upper, Nr);
        [m_NEW_LCR, ~] = LCR(STFT, IF_E, IM_E, sigma_s, cas);

        %% SNR signal
        X_win = 2*Lg:(L-2*Lg);
        for p = 1:Nr
            ref = modes(p, X_win);
            x_NEW_LCR = snr(ref, m_NEW_LCR(p, X_win) - ref);
            x_NEW_s = snr(ref, m_NEW_simple(p, X_win) - ref);
            x_S_RD = snr(ref, m_C_RD(p, X_win) - ref);
            x_MB_RD = snr(ref, m_MB_RD(p, X_win) - ref);
            if isnan(x_NEW_LCR + x_NEW_s + x_S_RD + x_MB_RD)
                error('MR : one of the SNR is NaN');
            end
        
            SNR_NEW_LCR(p, k) = SNR_NEW_LCR(p, k) + x_NEW_LCR;
            SNR_NEW_s(p, k) = SNR_NEW_s(p, k) + x_NEW_s;
            SNR_S_RD(p, k) = SNR_S_RD(p, k) + x_S_RD;
            SNR_MB_RD(p, k) = SNR_MB_RD(p, k) + x_MB_RD;
        end
        
        %% SNR IF
        for p = 1:Nr
            ref = IFs(p, X_win);
            x_NEW = snr(ref, IF_E(p, X_win) - ref);
            x_S_RD = snr(ref, L/Nfft*(Cs_simple(p, X_win) -1) - ref);
            x_MB_RD = snr(ref, L/Nfft*(Cs_VFB_MB(p, X_win) -1) - ref);
            if isnan(x_NEW + x_S_RD + x_MB_RD)
                error('IF : one of the SNR is NaN');
            end
        
            SNR_IF_NEW(p, k) = SNR_IF_NEW(p, k) + x_NEW;
            SNR_IF_S_RD(p, k) = SNR_IF_S_RD(p, k) + x_S_RD;
            SNR_IF_MB_RD(p, k) = SNR_IF_MB_RD(p, k) + x_MB_RD;
        end
    end
end

SNR_NEW_LCR = SNR_NEW_LCR/NRep;
SNR_NEW_s = SNR_NEW_s/NRep;
SNR_S_RD = SNR_S_RD/NRep;
SNR_MB_RD = SNR_MB_RD/NRep;

SNR_IF_NEW = SNR_IF_NEW/NRep;
SNR_IF_S_RD = SNR_IF_S_RD/NRep;
SNR_IF_MB_RD = SNR_IF_MB_RD/NRep;

end

