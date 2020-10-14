function [counts] = R1_TF_in_STD(s_in, L, Nfft, sigma_s, y, snr_in, k_IF_min, k_IF_max)

counts = zeros(3, 1);
noise = randn(L,1) + 1i*randn(L,1);
sn_LC = sigmerge(s_in, noise, snr_in);
[g, Lh] = create_gaussian_window(L, Nfft, sigma_s);

[STFT, omega, ~, QM, ~, tau] = FM_operators(sn_LC, Nfft, g, Lh, sigma_s);

[STFT_LM] = LM_from_STFT(STFT);

g_Vg = median(abs(real(STFT(:))))/0.6745;

Cg = y*g_Vg;

STFT_g = STFT.*(abs(STFT) > Cg);
STFT_LM_g = STFT_LM.*(abs(STFT) > Cg);

[id_RP_TFR, Energy_RP, E_RP_TFR, RP_maps] = R1_a_idRRP(STFT_LM_g, STFT, QM);
NB = length(Energy_RP);

[~, ~, ~, E_RB_TFR] = R1_b_idBasin(STFT_g, tau, omega, STFT_LM, id_RP_TFR, RP_maps, NB);

occ_k_vec = zeros(3, L);
for n=1:L
    [~, k_LM] = max(abs(STFT_LM_g(:, n)));
    [~, k_RRP] = max(E_RP_TFR(:, n));
    [~, k_RB] = max(E_RB_TFR(:, n));

    m_vec = [k_LM, k_RRP, k_RB];
    for p=1:3
        occ_k_vec(p, n) = m_vec(p);
        if m_vec(p) >= k_IF_min(n) && m_vec(p) <= k_IF_max(n)
            counts(p) = counts(p) + 1;
        end
    end
end

%% Figures
% figure;
% imagesc(abs(STFT_LM_g));
% set(gca,'ydir','normal');
% colormap(flipud(gray));
% title("STFT LM");
% hold on;
% plot(k_IF_min, 'r');
% plot(k_IF_max, 'r');
% plot(occ_k_vec(1, :));
% hold off;
% 
% figure;
% imagesc(E_RP_TFR);
% set(gca,'ydir','normal');
% colormap(flipud(gray));
% title("E RP");
% hold on;
% plot(occ_k_vec(2, :));
% hold off;
% 
% figure;
% imagesc(E_RB_TFR);
% set(gca,'ydir','normal');
% colormap(flipud(gray));
% title("E Basins");
% hold on;
% plot(occ_k_vec(3, :));
% hold off;

end

