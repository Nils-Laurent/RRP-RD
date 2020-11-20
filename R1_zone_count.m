function [counts] = R1_zone_count(STFT, QM, omega, tau, L, Nfft, y, k_IF_min, k_IF_max)

counts = zeros(3, 1);

[STFT_LM] = LM_from_STFT(STFT);

gamma_Vg = median(abs(real(STFT(:))))/0.6745;
C2_gamma = 2*gamma_Vg;

STFT_g = STFT.*(abs(STFT) > C2_gamma);
STFT_LM_g = STFT_LM.*(abs(STFT) > C2_gamma);

[id_RP_TFR, Energy_RP, E_RP_TFR, RP_maps] = R1_a_idRRP(STFT_LM_g, STFT, QM);
NB = length(Energy_RP);

[id_Basin_TFR, Energy_basin, ~, EB_RP_TFR] = R1_b_idBasin(STFT_g, tau, omega, STFT_LM, id_RP_TFR, RP_maps, NB);

[~, ~, ~, EZ_RP_TFR, ~] = R1_c_idZones(STFT_g, id_Basin_TFR, Energy_basin, id_RP_TFR);
%ref_Zone_TFR = EZ_RP_TFR.*(abs(STFT_g) > 0);
ref_Energy_TFR = EZ_RP_TFR;

for n=1:L
    [~, k_LM] = max(abs(STFT_LM_g(:, n)));
    
    k1 = k_IF_min(n);
    k2 = k_IF_max(n);
    
    Z_n = ref_Energy_TFR(:, n);
    k_vec_zone = Z_n == max(Z_n);
    if sum(k_vec_zone(k1:k2)) > 0
        counts(3) = counts(3) + 1;
    end
    
    if k_LM >= k1 && k_LM <= k2
        counts(2) = counts(2) + 1;
        counts(1) = counts(1) + 1;
        continue;
    end
    c_n = sum(abs(STFT_LM_g(k1:k2, n)));
    if c_n > 0
        counts(1) = counts(1) + 1;
    end
end

%% Figures
% figure;
% imagesc(abs(STFT_LM_g));
% set(gca,'ydir','normal');
% colormap(flipud(gray));
% title("STFT LM");
% hold on;
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
% imagesc(EB_RP_TFR);
% set(gca,'ydir','normal');
% colormap(flipud(gray));
% title("E Basins RP");
% hold on;
% plot(occ_k_vec(3, :));
% hold off;

% --- log
% figure;
% imagesc(log(abs(STFT_LM_g)));
% set(gca,'ydir','normal');
% colormap(flipud(gray));
% title("STFT LM (log)");
% hold on;
% plot(occ_k_vec(1, :));
% hold off;
% 
% figure;
% imagesc(log(E_RP_TFR));
% set(gca,'ydir','normal');
% colormap(flipud(gray));
% title("E RP (log)");
% hold on;
% plot(occ_k_vec(2, :));
% hold off;
% 
% figure;
% imagesc(log(EB_RP_TFR));
% set(gca,'ydir','normal');
% colormap(flipud(gray));
% title("E Basins RP (log)");
% hold on;
% plot(occ_k_vec(3, :));
% hold off;

end

