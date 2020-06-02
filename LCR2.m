function [modes, TFR_denoised, Lg] = LCR2(STFT, g, Lg, Nfft, XCs, IFs, IMs, sigma_s)
[~, L] = size(STFT);
[Nr, ~] = size(IFs);

% [g, Lg] = create_gaussian_window(L, Nfft, sigma_s);


%% MR
TFR_denoised = zeros(size(STFT));
modes = zeros(Nr, length(XCs));

for p = 1:Nr
    %% use estimate and inverse STFT
     [TFR_denoised_r] = LCR_estim_STFT(sigma_s, STFT, IFs(p, :), IMs(p, :), Nfft, XCs);
    TFR_denoised = TFR_denoised + TFR_denoised_r;
%     modes(p, :) = L*itfrstft(TFR_denoised_r(:, XCs), cas, g);
    
    %case without periodizing
    x = zeros(1,length(XCs));
    Lg = (length(g)-1)/2;
    for n=1:length(XCs)
        icol = XCs(n);
        x(n) = L/g(Lg+1)*sum(TFR_denoised_r(:,icol))/Nfft;
    end
    modes(p, :) = x;
end

% figure;
% imagesc(abs(TFR_denoised));
% % imagesc((0:L-1)*T/L, (0:Nfft-1)*L/(Nfft*T), abs(TFR));
% set(gca,'ydir','normal');
% colormap(flipud(gray));
% axis square
% pause;

end
