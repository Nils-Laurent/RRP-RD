function [TFR_denoised] = LCR_estim_STFT(sigma_s, STFT, phipE, phippE, Nfft, XCs)

[N_STFT, L] = size(STFT);
TFR_denoised = zeros(size(STFT));
for nCs = 1:length(XCs)
    n = XCs(nCs);
    Zn = pi*sigma_s^2*(1 + 1i*phippE(nCs)*sigma_s^2)...
        /(1 + phippE(nCs)^2*sigma_s^4);
    
    KWd = round(phipE(nCs)*Nfft/L)+1;
    KWd = min(N_STFT, KWd);
    KWd = max(1, KWd);
    
    Wd = ((KWd-1)*L/Nfft);
    Gn = STFT(KWd, n)*exp(Zn*(Wd - phipE(nCs))^2);
    
    %% assign coefficients
    for k = 1:N_STFT
        TFR_denoised(k, n) = Gn*exp(-Zn*((k-1)*L/Nfft - phipE(nCs))^2);
    end
end

end

