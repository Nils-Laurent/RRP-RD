function [E_P_vec] = R1_prop_E(STFT, QM, omega, tau, L, Nfft, P_max)

[STFT_LM] = LM_from_STFT(STFT);

gamma_Vg = median(abs(real(STFT(:))))/0.6745;

C2_gamma = 2*gamma_Vg;

STFT_g = STFT.*(abs(STFT) > C2_gamma);
STFT_LM_g = STFT_LM.*(abs(STFT) > C2_gamma);

[id_RP_TFR, Energy_RP, ~, RP_maps] = R1_a_idRRP(STFT_LM_g, STFT, QM);
NB = length(Energy_RP);

[id_Basin_TFR, Energy_basin, ~, ~] = R1_b_idBasin(STFT_g, tau, omega, STFT_LM, id_RP_TFR, RP_maps, NB);

[id_Zone_TFR, Energy_zone, idZ_RP_TFR, EZ_RP_TFR, Zone_E_TFR] = R1_c_idZones(STFT_g, id_Basin_TFR, Energy_basin, id_RP_TFR);
NZ = length(Energy_zone);

[~, E_P_vec] = R1_d_selection(idZ_RP_TFR, EZ_RP_TFR, NB, NZ, P_max);
end

