function [] = R1_RRP_RD(STFT, QM, omega, tau, L, Nfft, P, sigma_s)

[STFT_LM] = LM_from_STFT(STFT);

gamma_Vg = median(abs(real(STFT(:))))/0.6745;

C2_gamma = 2*gamma_Vg;
C3_gamma = 3*gamma_Vg;

ASTFT_g2 = abs(STFT).*(abs(STFT) > C2_gamma);
A_LM_g2 = abs(STFT_LM).*(abs(STFT) > C2_gamma);
A_LM_g3 = abs(STFT_LM).*(abs(STFT) > C3_gamma);

[id_RP_TFR, Energy_RP, E_RP_TFR, RP_maps] =...
    R1_a_idRRP(A_LM_g2, STFT, QM);
NB = length(Energy_RP);

[id_Basins_TFR, Energy_basins, E2_basins, E_Basins_TFR, E2_Basins_TFR, EB_RP_TFR] =...
    R1_b_idBasin(ASTFT_g2, A_LM_g3, tau, omega, id_RP_TFR, RP_maps, NB);

[id_Zones_TFR, Energy_Zones, idZ_RP_TFR, EZ_RP_TFR, E_Zones_TFR] =...
    R1_c_idZones(ASTFT_g2, id_Basins_TFR, Energy_basins, id_RP_TFR);

%% new data based on g3
EZ_new = zeros(size(Energy_Zones));
EB_new = zeros(size(Energy_basins));
for n=1:L
    for k=1:Nfft
        e_g3 = A_LM_g3(k, n);
        if e_g3 > 0
            id_zone = idZ_RP_TFR(k, n);
            EZ_new(id_zone) = EZ_new(id_zone) + e_g3;
            
            id_basin = id_Basins_TFR(k, n);
            if id_basin > 0
                EB_new(id_basin) = EB_new(id_basin) + e_g3;
            end
        end
    end
end

EZ_new_TFR = zeros(size(E_Zones_TFR));
EZ_new_RP_TFR = zeros(size(E_Zones_TFR));
EB_new_TFR = zeros(size(E_Basins_TFR));
for n=1:L
    for k=1:Nfft
        id_zone = id_Zones_TFR(k, n);
        if id_zone > 0
            EZ_new_TFR(k, n) = EZ_new(id_zone);
        end
        
        id_zone_RP = idZ_RP_TFR(k, n);
        if id_zone_RP > 0
            EZ_new_RP_TFR(k, n) = EZ_new(id_zone_RP);
        end
        
        id_basin = id_Basins_TFR(k, n);
        if id_basin > 0
            EB_new_TFR(k, n) = EB_new(id_basin);
        end
    end
end
E2_basins = EB_new;
E2_Basins_TFR = EB_new_TFR;
Energy_Zones = EZ_new;

%% figures
% figure;
% imagesc(abs(A_LM_g2));
% set(gca,'ydir','normal');
% colormap(flipud(gray));
% title("LM thresh (2)");
% pbaspect([1 1 1]);
% set(gcf, 'Position',  [0, 0, 1000, 1000]);
% 
% figure;
% imagesc(abs(A_LM_g3));
% set(gca,'ydir','normal');
% colormap(flipud(gray));
% title("LM thresh (3)");
% pbaspect([1 1 1]);
% set(gcf, 'Position',  [0, 0, 1000, 1000]);

% max_vec = zeros(2, L);
% for n=1:L
%     E_vec = E_Zones_TFR(:, n);
%     [ES_vec, s_pos] = sort(E_vec, 'descend');
%     [~, u_pos] = unique(nonzeros(ES_vec), 'stable');
%     if length(u_pos) >= 1
%         max_vec(1, n) = s_pos(u_pos(1));
%     end
%     if length(u_pos) >= 2
%         max_vec(2, n) = s_pos(u_pos(2));
%     end
% end
% 
% figure;
% imagesc(E_Zones_TFR);
% set(gca,'ydir','normal');
% colormap(flipud(gray));
% title("Energy zones");
% hold on;
% plot(max_vec(1, :));
% plot(max_vec(2, :));
% hold off;
% pbaspect([1 1 1]);
% set(gcf, 'Position',  [0, 0, 1000, 1000]);

% figure;
% imagesc(E_Basins_TFR);
% set(gca,'ydir','normal');
% colormap(flipud(gray));
% title("Energy basins");
% pbaspect([1 1 1]);
% set(gcf, 'Position',  [0, 0, 1000, 1000]);

% figure;
% imagesc(E2_Basins_TFR);
% set(gca,'ydir','normal');
% colormap(flipud(gray));
% title("Energy2 basins");
% pbaspect([1 1 1]);
% set(gcf, 'Position',  [0, 0, 1000, 1000]);

% return;
% return;
% close all;
R1_e_spline(id_Basins_TFR, id_Zones_TFR, idZ_RP_TFR, E2_Basins_TFR, A_LM_g3, E2_basins, Energy_Zones, P);
% R1_e_spline(id_Zones_TFR, idZ_RP_TFR, EZ_new_TFR, EZ_new_RP_TFR, EB_new_TFR, STFT, A_LM_g3, EZ_new, NB, sigma_s, P);

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

