function [Modes_max, E_max] = R1_RRP_RD(STFT, QM, omega, tau, Fs, Nfft, Nr, sigma_s, smooth_p)

% STFT : Short time fourier transform
% QM : 2nd order modulation operator
% omega : IF estimate
% tau : groupe delay
% Fs : Sampling frequency
% Nfft : Number of frequency bins
% Nr : Number of modes
% sigma_s (unused) : gaussian coefficient
% smooth_p : smoothing quantity to apply on the spline approximation

[STFT_LM] = LM_from_STFT(STFT);

gamma_Vg = median(abs(real(STFT(:))))/0.6745;

C2_gamma = 2*gamma_Vg;
C3_gamma = 3*gamma_Vg;

ASTFT_g2 = abs(STFT).*(abs(STFT) > C2_gamma);
A_LM_g2 = abs(STFT_LM).*(abs(STFT) > C2_gamma);
A_LM_g3 = abs(STFT_LM).*(abs(STFT) > C3_gamma);

[id_RP_TFR, Energy_RP, E_RP_TFR, RP_maps] =...
    R1_a_idRRP(A_LM_g2, STFT, QM, Nfft, Fs);
NB = length(Energy_RP);

% figure;
% imagesc(A_LM_g2);
% set(gca,'ydir','normal');
% colormap(flipud(gray));
% title("LM g2");
% 
% figure;
% imagesc(A_LM_g3);
% set(gca,'ydir','normal');
% colormap(flipud(gray));
% title("LM g3");
% pause;

[id_Basins_TFR, Energy_basins, E2_basins, E_Basins_TFR, E2_Basins_TFR, EB_RP_TFR] =...
    R1_b_idBasin(ASTFT_g2, A_LM_g3, tau, omega, id_RP_TFR, RP_maps, NB, Nfft, Fs);

% figure;
% imagesc(E_Basins_TFR);
% set(gca,'ydir','normal');
% colormap(flipud(gray));
% title("E basins");

[id_Zones_TFR, Energy_Zones, idZ_RP_TFR, EZ_RP_TFR, E_Zones_TFR] =...
    R1_c_idZones(ASTFT_g2, id_Basins_TFR, Energy_basins, id_RP_TFR);

% figure;
% imagesc(E_Zones_TFR);
% set(gca,'ydir','normal');
% colormap(flipud(gray));
% title("E zones");
% pause;

% fprintf("Set splines\n");

[Modes_max, E_max] = R1_e_spline(E_Zones_TFR, id_Basins_TFR, id_Zones_TFR,...
    idZ_RP_TFR, E2_Basins_TFR, A_LM_g3, E2_basins, Energy_Zones, Nr, smooth_p, Fs, Nfft);

return;

%% analyse des zones

% g_basin_disp = id_Basin_TFR + 14*(id_Basin_TFR > 0);
% figure;
% imagesc(g_basin_disp);
% set(gca,'ydir','normal');
% % myColorMap = colorcube(round(NB/8));
% myColorMap = colorcube;
% myColorMap(1,:) = 1;
% colormap(myColorMap);
% title("Basin ID");
% 
% g_zone_disp = id_Zone_TFR.*(id_Zone_TFR > NB);
% g_zone_disp = g_zone_disp - NB*(id_Zone_TFR > NB);
% 
% figure;
% imagesc(g_zone_disp);
% set(gca,'ydir','normal');
% myColorMap = lines(NZ - NB);
% myColorMap(1,:) = 1;
% colormap(myColorMap);
% title("Zone ID");
% return;

end

