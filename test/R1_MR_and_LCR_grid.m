function [m_SR, m_LCR, IF_vecs, STFT_LCR] = R1_MR_and_LCR_grid(STFT, QM, Cs, g, Lh, sigma_s, Nr, Nfft, Fs)

[~, L] = size(STFT);

IF_vecs = zeros(Nr, L);

t = (0:L-1)/Fs;
IM_vecs = zeros(Nr, L);
K_lower = zeros(Nr, L);
K_upper = zeros(Nr, L);
Range_eta = zeros(Nr, L);

for p=1:Nr
%     IF_vecs(p, :) = fnval(Spl(p).spline, t);
    IF_vecs(p, :) = (Cs(p, :) - 1)*Fs/Nfft;
%     IM_vecs(p, :) = fnval(fnder(Spl(p).spline), t);
    for n=1:L
        IM_vecs(p, n) = real(QM(Cs(p, n), n));
    end
    K_mid = round(IF_vecs(p, :)*Nfft/Fs) + 1;
    Range_eta(p, :) = 1/(sqrt(2*pi)*sigma_s)*sqrt(1 + sigma_s^4*IM_vecs(p, :).^2);
    range_p = ceil(3*Range_eta(p, :)*Nfft/Fs);
    K_lower(p, :) = K_mid - range_p;
    K_upper(p, :) = K_mid + range_p;
end

% [N_Y, ~] = size(STFT);
% figure;
% imagesc(t, (0:N_Y-1)*Fs/Nfft, abs(STFT));
% set(gca,'ydir','normal');
% colormap(flipud(gray));
% axis square;
% hold on;
% for m_spl=1:Nr
%     plot(t, IF_vecs(m_spl, :) - Range_eta(p, :), '--');
%     plot(t, IF_vecs(m_spl, :));
%     plot(t, IF_vecs(m_spl, :) + Range_eta(p, :), '--');
% end
% hold off;
% pause;

[m_SR, ~] = MR_simple(STFT, Fs, Nfft, g, Lh, K_lower, K_upper, Nr);

cas = 1;
[m_LCR, STFT_LCR] = LCR(STFT, IF_vecs, IM_vecs, sigma_s, Fs, Nfft, cas);

end

