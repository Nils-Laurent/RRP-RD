function [modes, TFR_denoised, Lg] = LCR(STFT, IFs, IMs, sigma_s, Fs, Nfft, cas)
[N_Y, L] = size(STFT);
[Nr, ~] = size(IFs);

[g, Lg] = create_gaussian_window(L, Nfft, sigma_s);


%% MR
TFR_denoised = zeros(size(STFT));
modes = zeros(Nr, L);

for p = 1:Nr
    %% use estimate and inverse STFT
     [TFR_denoised_r] = LCR_estim_STFT(sigma_s, STFT, IFs(p, :), IMs(p, :), Nfft, Fs);
    TFR_denoised = TFR_denoised + TFR_denoised_r;
    % modes(p, :) = L*itfrstft(TFR_denoised_r, cas, g);
    modes(p, :) = FM_inverse(TFR_denoised_r, Fs, Nfft, g, cas);
end

end
