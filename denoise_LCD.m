function [s_denoised_LC, TFR_denoised, Cs, Lg, E2] = denoise_LCD(s_noise, NRidges, sigma_s, Nfft)

L = length(s_noise);
cas = 1;

[g, Lg] = create_gaussian_window(L, Nfft, sigma_s);


[TFR_noise] = tfrstft(s_noise, Nfft, cas, g, Lg);
TFR_noise = TFR_noise/L;

%% ridge extraction
[Cs_1] = exridge_mult(TFR_noise, NRidges, 0, 0, 2);
[TFR_n, omega, omega2, q] = FM_operators(s_noise, Nfft, g, Lg, sigma_s);
[Cs] = exridge_new(TFR_n, Lg, sigma_s, q, omega, omega2, 2);

for n=1:L
    if isnan(Cs(n))
        Cs(n) = Cs_1(n);
    end
end
Cs = Cs';
% Cs = exridge(TFR_noise, 0, 0, 2);

TFR_denoised = zeros(size(TFR_noise));
E2 = zeros(NRidges, 2, L);
%gamma = median(abs(real(TFR_noise(:))))/0.6745;
for r = 1:NRidges
    %% mode estimation
    [phipE1, phipE2, phippE] = retrieve_mode(s_noise, Nfft, g, Lg, sigma_s, Cs(r, :));
    X = [phipE1, phipE2];
    E2(r, :, :) = transpose(X);

    %% use estimate and inverse STFT
    [TFR_denoised_r] = tfr_from_estimation(sigma_s, TFR_noise, phipE2, phippE, L, Nfft);
    TFR_denoised = TFR_denoised + TFR_denoised_r;
end

[s_denoised_LC] = L*itfrstft(TFR_denoised, cas, g);

% figure;
% imagesc((0:L-1)/L, (L/Nfft)*(1:Nfft), abs(TFR_denoised));
% set(gca,'ydir','normal');
% axis square
% colorbar;
% pause;

% XLg = 2*Lg:L-2*Lg;
% TFR_D_mod = abs(TFR_noise(:, XLg)) - abs(TFR_denoised(:, XLg));
% TFR_D_real = real(TFR_noise(:, XLg)) - real(TFR_denoised(:, XLg));
% TFR_D_imag = imag(TFR_noise(:, XLg)) - imag(TFR_denoised(:, XLg));
% 
% TFR_D_mod_all = abs(TFR_noise) - abs(TFR_denoised);
% 
% figure;
% imagesc(abs(TFR_D_real) - abs(TFR_D_imag));
% title("|real error| - |imag error|");
% set(gca,'ydir','normal');
% colorbar;
% figure;
% imagesc(XLg, 1:Nfft, TFR_D_mod);
% title("modulus error");
% set(gca,'ydir','normal');
% colorbar;

% XERR = zeros(L, 1);
% for index = 1:L
%     kRidge = round(phip(index)*Nfft/L)+1;
%     XERR(index) = TFR_D_mod_all(kRidge, index);
% end
% 
% figure;
% F_XERR = fft(XERR, Nfft);
% F_XERR = F_XERR(Nfft/2+1:Nfft);
% plot((0:Nfft/2-1)*(L/Nfft), abs(F_XERR));
% title("error frequency");
% [fval, freq] = max(F_XERR);
% freq = (freq-1)*L/Nfft;
% fprintf("max|fft(Xerr)| at %d Hz\n", freq);
% factor(freq)
% 
% T0 = 1770;
% figure;
% subplot(2, 1, 1);
% plot(1:Nfft,real(TFR_noise(:,T0)),'k', 1:Nfft,real(TFR_denoised(:,T0)),'b--');
% title(sprintf("real, T0 = %d", T0));
% subplot(2, 1, 2);
% plot(1:Nfft,imag(TFR_noise(:,T0)),'k', 1:Nfft,imag(TFR_denoised(:,T0)),'r--');
% title(sprintf("imag, T0 = %d", T0));
% figure;
% plot(1:Nfft,real(TFR_denoised(:,T0)) - real(TFR_noise(:,T0)), 1:Nfft,imag(TFR_denoised(:,T0)) - imag(TFR_noise(:,T0)), '--');
% title(sprintf("errors, T0 = %d", T0));
% % pause

end
