function [modes, TFR_denoised, Lg] = LCR_partial(STFT, g, Lg, Nfft, XCs, IFs, IMs, sigma_s)
[~, L] = size(STFT);
[Nr, ~] = size(IFs);

%% MR
TFR_denoised = zeros(size(STFT));
modes = zeros(Nr, length(XCs));

for p = 1:Nr
    %% use estimate and inverse STFT
     [TFR_denoised_r] = LCR_estim_STFT(sigma_s, STFT, IFs(p, :), IMs(p, :), Nfft, XCs);
    TFR_denoised = TFR_denoised + TFR_denoised_r;
    
    %case without periodizing
    x = zeros(1,length(XCs));
    Lg = (length(g)-1)/2;
    for n=1:length(XCs)
        icol = XCs(n);
        x(n) = L/g(Lg+1)*sum(TFR_denoised_r(:,icol))/Nfft;
    end
    modes(p, :) = x;
end

end
