function [RE_vec] = sig_min(sigma_set, s_clean, L, Nfft)
    SL = length(sigma_set);
    RE_vec = zeros(1, SL);
    iSL = 0;
    alpha = 3;
    
    for sigma = sigma_set
        iSL = iSL + 1;
        fprintf('%u/%u\n', iSL, SL);
        [g, Lg] = create_gaussian_window(L, Nfft, sigma);
        [TFR, ~, ~, ~] = FM_operators(s_clean, Nfft, g, Lg, sigma);
        Y = abs(TFR);
        TFR_MS = sum(Y(:));
        RE_vec(iSL) = 1/(1 - alpha)*log2(sum(Y(:).^alpha)/(TFR_MS^alpha)) - log2(L);
    end
    
end