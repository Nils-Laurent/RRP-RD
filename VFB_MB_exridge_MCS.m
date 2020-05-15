function [Cs] = VFB_MB_exridge_MCS(STFT, sigma_s, QM, C, Nr)
    % C = 2
    % -round(Nfft*Q/N^2)
    
    [Nfft, L] = size(STFT);
    
    Cs = zeros(Nr, L);
    
    RK = -round(Nfft*real(QM)/L^2);
    STFT_it = STFT;
    for p=1:Nr
        [Cr, ~] = VFB_MB_exridge(STFT_it, 0, 0, RK, C);
        Cs(p, :) = Cr;
        [~, STFT_it] = VFB_MB_extrae_modo(STFT_it, Cr, real(QM), L, Nfft, sigma_s);
    end
end

