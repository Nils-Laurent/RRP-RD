function [modes, TFR_denoised, Lg] = LCR(STFT, IFs, IMs, sigma_s, cas)
[Nfft, L] = size(STFT);
[Nr, ~] = size(IFs);

[g, Lg] = create_gaussian_window(L, Nfft, sigma_s);


%% MR
TFR_denoised = zeros(size(STFT));
modes = zeros(Nr, L);

for p = 1:Nr
    %% use estimate and inverse STFT
     [TFR_denoised_r] = LCR_estim_STFT(sigma_s, STFT, IFs(p, :), IMs(p, :), Nfft, 1:L);
    TFR_denoised = TFR_denoised + TFR_denoised_r;
    modes(p, :) = L*itfrstft(TFR_denoised_r, cas, g);
end

% figure;
% imagesc(abs(TFR_denoised));
% % imagesc((0:L-1)*T/L, (0:Nfft-1)*L/(Nfft*T), abs(TFR));
% set(gca,'ydir','normal');
% colormap(flipud(gray));
% axis square
% pause;

end
