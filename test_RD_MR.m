function [SNR_NEW, SNR_C_RD, SNR_MB_RD,...
    SNR_IF_NEW, SNR_IF_C_RD, SNR_IF_MB_RD] = test_RD_MR(...
    modes, IFs, clwin, sigma_s, Nfft, poly_d, SNR_IN, NRep)

[Nr, L] = size(modes);
s_in = sum(modes, 1);

SNR_NEW = zeros(Nr, length(SNR_IN));
SNR_C_RD = zeros(Nr, length(SNR_IN));
SNR_MB_RD = zeros(Nr, length(SNR_IN));

SNR_IF_NEW = zeros(Nr, length(SNR_IN));
SNR_IF_C_RD = zeros(Nr, length(SNR_IN));
SNR_IF_MB_RD = zeros(Nr, length(SNR_IN));

for k=1:length(SNR_IN)
    for l=1:NRep
        fprintf('snr %d/%d, rep %d/%d\n', k, length(SNR_IN), l, NRep);
        
        noise = randn(L, 1)+1i*randn(L, 1);
        s_noise = sigmerge(transpose(s_in), noise, SNR_IN(k));
        [g, Lg] = create_gaussian_window(L, Nfft, sigma_s);
        [STFT, ~, ~, QM] = FM_operators(s_noise, Nfft, g, Lg, sigma_s);

%         figure;
%         imagesc(1:L, 1:Nfft, abs(STFT));
%         set(gca,'ydir','normal');
%         axis square
%         colormap(flipud(gray));
%         title("STFT");
%         pause;
        
        fprintf('VFB MB, ');
        [Cs_VFB_MB] = VFB_MB_exridge_MCS(STFT, sigma_s, QM, 2, Nr);
        [K_lower, K_upper] = MR_interval(Cs_VFB_MB, QM, Nfft, sigma_s);
        [m_MB_RD, ~] = MR_simple(STFT, g, Lg, K_lower, K_upper, Nr);
        
        fprintf('classic, ');
        [Cs_classic] = exridge_mult(STFT, Nr, 0, 0, clwin);
        [K_lower, K_upper] = MR_interval(Cs_classic, QM, Nfft, sigma_s);
        [m_C_RD, ~] = MR_simple(STFT, g, Lg, K_lower, K_upper, Nr);
        
        fprintf('RRP RD\n');
        [~, Qs, KY_lower, KY_upper, ~] = novel_RRP_RD(STFT, QM, sigma_s, Nr, poly_d);
        Cs = zeros(size(modes));
        for p = 1:Nr
            Cs(p, :) = round(polyval(Qs(p, :), (0:L-1)/L)*Nfft/L) + 1;
            Cs(p, :) = min(Nfft, max(1, Cs(p, :)));
        end
        [m_NEW, ~] = MR_simple(STFT, g, Lg, KY_lower, KY_upper, Nr);

        
        %% SNR signal
        X_win = 2*Lg:(L-2*Lg);
        for p = 1:Nr
            ref = modes(p, X_win);
            x_NEW = snr(ref, m_NEW(p, X_win) - ref);
            x_C_RD = snr(ref, m_C_RD(p, X_win) - ref);
            x_MB_RD = snr(ref, m_MB_RD(p, X_win) - ref);
        
            SNR_NEW(p, k) = SNR_NEW(p, k) + x_NEW;
            SNR_C_RD(p, k) = SNR_C_RD(p, k) + x_C_RD;
            SNR_MB_RD(p, k) = SNR_MB_RD(p, k) + x_MB_RD;
        end
        
        %% SNR IF
        for p = 1:Nr
            ref = IFs(p, X_win);
            x_NEW = snr(ref, Cs(p, X_win) - ref);
            x_C_RD = snr(ref, Cs_classic(p, X_win) - ref);
            x_MB_RD = snr(ref, Cs_VFB_MB(p, X_win) - ref);
        
            SNR_IF_NEW(p, k) = SNR_IF_NEW(p, k) + x_NEW;
            SNR_IF_C_RD(p, k) = SNR_IF_C_RD(p, k) + x_C_RD;
            SNR_IF_MB_RD(p, k) = SNR_IF_MB_RD(p, k) + x_MB_RD;
        end
    end
end

SNR_NEW = SNR_NEW/NRep;
SNR_C_RD = SNR_C_RD/NRep;
SNR_MB_RD = SNR_MB_RD/NRep;

SNR_IF_NEW = SNR_IF_NEW/NRep;
SNR_IF_C_RD = SNR_IF_C_RD/NRep;
SNR_IF_MB_RD = SNR_IF_MB_RD/NRep;

end

