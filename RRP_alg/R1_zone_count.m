function [counts] = R1_zone_count(STFT, QM, omega, tau, L, Nfft, y, k_IF_min, k_IF_max)

counts = zeros(3, 1);

[STFT_LM] = LM_from_STFT(STFT);

gamma_Vg = median(abs(real(STFT(:))))/0.6745;

C2_gamma = 2*gamma_Vg;
C3_gamma = 3*gamma_Vg;

ASTFT_g2 = abs(STFT).*(abs(STFT) > C2_gamma);
A_LM_g2 = abs(STFT_LM).*(abs(STFT) > C2_gamma);
A_LM_g3 = abs(STFT_LM).*(abs(STFT) > C3_gamma);

[id_RP_TFR, Energy_RP, E_RP_TFR, RP_maps] =...
    R1_a_idRRP(A_LM_g2, STFT, QM, Nfft, L);
NB = length(Energy_RP);

[id_Basins_TFR, Energy_basins, E2_basins, E_Basins_TFR, E2_Basins_TFR, EB_RP_TFR] =...
    R1_b_idBasin(ASTFT_g2, A_LM_g3, tau, omega, id_RP_TFR, RP_maps, NB, Nfft, L);

[id_Zones_TFR, Energy_Zones, idZ_RP_TFR, EZ_RP_TFR, E_Zones_TFR] =...
    R1_c_idZones(ASTFT_g2, id_Basins_TFR, Energy_basins, id_RP_TFR);

ref_Energy_TFR = EZ_RP_TFR;

for n=1:L
    [~, k_LM] = max(abs(A_LM_g2(:, n)));
    
    k1 = k_IF_min(n);
    k2 = k_IF_max(n);
    
    Z_n = ref_Energy_TFR(:, n);
    k_vec_zone = Z_n == max(Z_n);
    if sum(k_vec_zone(k1:k2)) > 0
        % Zone global max : Q(beta)
        counts(3) = counts(3) + 1;
    end
    
    if k_LM >= k1 && k_LM <= k2
        % LM global max and existence of LM
        counts(2) = counts(2) + 1;
        counts(1) = counts(1) + 1;
        continue;
    end
    c_n = sum(abs(A_LM_g2(k1:k2, n)));
    if c_n > 0
        % existence of LM : P(beta)
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

