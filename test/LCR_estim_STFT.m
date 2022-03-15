function [TFR_denoised] = LCR_estim_STFT(sigma_s, STFT, phipE, phippE, Nfft, Fs)

[N_STFT, L] = size(STFT);
TFR_denoised = zeros(size(STFT));
for n = 1:L
    Zn = pi*sigma_s^2*(1 + 1i*phippE(n)*sigma_s^2)...
        /(1 + phippE(n)^2*sigma_s^4);
    
    KWd = round(phipE(n)*Nfft/Fs)+1;
    KWd = min(N_STFT, KWd);
    KWd = max(1, KWd);
    
    Wd = ((KWd-1)*Fs/Nfft);
    Gn = STFT(KWd, n)*exp(Zn*(Wd - phipE(n))^2);
    
    %% assign coefficients
    for k = 1:N_STFT
        TFR_denoised(k, n) = Gn*exp(-Zn*((k-1)*Fs/Nfft - phipE(n))^2);
    end
end

end

