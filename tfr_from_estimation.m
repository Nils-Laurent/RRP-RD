function [TFR_denoised] = tfr_from_estimation(sigma, STFT, phipE, phippE, L, Nfft)

TFR_denoised = zeros(Nfft, L);
for n = 1:L
    Zn = pi*sigma^2*(1 + 1i*phippE(n)*sigma^2)...
        /(1 + phippE(n)^2*sigma^4);
    
    KWd = round(phipE(n)*Nfft/L)+1;
    KWd = min(Nfft, KWd);
    KWd = max(1, KWd);
    
    Wd = ((KWd-1)*L/Nfft);
    Gn = STFT(KWd, n)*exp(Zn*(Wd - phipE(n))^2);
    
    %% assign coefficients
    for k = 1:Nfft
        TFR_denoised(k, n) = Gn*exp(-Zn*((k-1)*L/Nfft - phipE(n))^2);
    end
end

end

