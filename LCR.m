function [modes, TFR_denoised, Lg] = LCR(STFT, IFs, IMs, sigma_s, cas)
[Nfft, L] = size(STFT);
[Nr, ~] = size(IFs);

[g, Lg] = create_gaussian_window(L, Nfft, sigma_s);


%% MR
TFR_denoised = zeros(size(STFT));
modes = zeros(Nr, L);

for p = 1:Nr
%     m = 0;
%     iPX = [];
%     iPY = [];
%     iX = zeros(1, L);
%     for n=1:L
%         X = sigma_s*(1 + sigma_s^4*IMs(p, n)^2)^(-1/4);
%         iX(n) = X;
%         if D_RRP(p, n) == 0
%             continue;
%         end
%         m = m + 1;
%         iPX(m) = (n-1)/L;
%         iPY(m) = abs(STFT(D_RRP(p, n), n))/X;
%     end
%     
% %     degree_WPF = 5;
% %     PW = polyfit(iPX, iPY, degree_WPF);
% %     AE = polyval(PW, (0:L-1)/L);
%     AE = pchip(iPX, iPY, (0:L-1)/L);
%     
%     figure;
%     
%     subplot(3, 1, 1);
%     plot(iPX, iX);
%     title("iX");
%     
%     subplot(3, 1, 2);
%     plot(iPX, iPY);
%     title("iPY");
%     
%     subplot(3, 1, 3);
%     plot(AE);
%     title("AE");
    

    %% use estimate and inverse STFT
    %[TFR_denoised_r] = tfr_from_estimation(sigma_s, STFT, phipE2, phippE, L, Nfft);
     [TFR_denoised_r] = LCR_estim_STFT(sigma_s, STFT, IFs(p, :), IMs(p, :), L, Nfft);
    TFR_denoised = TFR_denoised + TFR_denoised_r;
    modes(p, :) = L*itfrstft(TFR_denoised_r, cas, g);
end

% figure;
% imagesc(1:L, 1:Nfft, abs(STFT));
% set(gca,'ydir','normal');
% axis square
% colormap(flipud(gray));
% title("STFT noise");
% 
% figure;
% imagesc(1:L, 1:Nfft, abs(TFR_denoised));
% set(gca,'ydir','normal');
% axis square
% colormap(flipud(gray));
% title("STFT LCR");

end
