function [E_P_vec] = R1_prop_E(STFT, QM, omega, tau, L, Nfft, P_max)

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
NZ = length(Energy_Zones);

[~, E_P_vec] = R1_d_selection(idZ_RP_TFR, EZ_RP_TFR, NB, NZ, P_max);

end

