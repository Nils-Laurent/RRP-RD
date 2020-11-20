function [id_Zone_TFR, Energy_zone, idZ_RP_TFR, EZ_RP_TFR, E_Zone_TFR] =...
    R1_c_idZones(STFT_th, id_Basin_TFR, Energy_basin, id_Ridge_TFR)

[Nfft, L] = size(STFT_th);

Basin_ID_TH_TFR = id_Basin_TFR.*(abs(STFT_th) > 0);
N_basin = length(Energy_basin);

map_basin = zeros(1, N_basin);
N_zone = N_basin;


% figure;
% imagesc(abs(Basin_ID_TH_TFR));
% set(gca,'ydir','normal');
% colormap(flipud(gray));
% axis square;
% title("ID to zone");
% pbaspect([1 1 1]);
% set(gcf, 'Position',  [0, 0, 1000, 1000]);

%% v4
id_vec = zeros(3, 1);
for n=1:(L - 1)
    for k=1:(Nfft - 1)
        a = Basin_ID_TH_TFR(k, n);
        if a == 0
            continue;
        end
        
        id_vec(1) = Basin_ID_TH_TFR(k+1, n);
        id_vec(2) = Basin_ID_TH_TFR(k, n+1);
        id_vec(3) = Basin_ID_TH_TFR(k+1, n+1);
        
        for p=1:3
            b = id_vec(p);
            if b == 0
                continue;
            end
            
            za = map_basin(a);
            zb = map_basin(b);
            z_min = min(za, zb);
            if z_min > 0
%                 z_id = max(za, zb);
                map_basin(map_basin == z_min) = max(za, zb);
            else
                z_id = max(za, zb);
                if z_id == 0
                    N_zone = N_zone + 1;
                    z_id = N_zone;
                end
                map_basin(a) = z_id;
                map_basin(b) = z_id;
            end
        end
    end
end

%% create output data
Energy_zone = zeros(N_zone, 1);
for basin_id=1:N_basin
    if map_basin(basin_id) == 0
        map_basin(basin_id) = basin_id;
    end
    zone_id = map_basin(basin_id);
    Energy_zone(zone_id) = Energy_zone(zone_id) + Energy_basin(basin_id);
end

id_Zone_TFR = zeros(size(STFT_th));
E_Zone_TFR = zeros(size(STFT_th));
idZ_RP_TFR = zeros(size(STFT_th));
EZ_RP_TFR = zeros(size(STFT_th));
for n=1:L
    for k=1:Nfft
        basin_id = id_Basin_TFR(k, n);
        if basin_id > 0
            zone_id = map_basin(basin_id);
            id_Zone_TFR(k, n) = zone_id;
            E_Zone_TFR(k, n) = Energy_zone(zone_id);
        end
        
        RP_id = id_Ridge_TFR(k, n);
        if RP_id > 0
            zone_id = map_basin(RP_id);
            idZ_RP_TFR(k, n) = zone_id;
            EZ_RP_TFR(k, n) = Energy_zone(map_basin(RP_id));
        end
    end
end

return;
%% figures A
figure;
imagesc(id_Zone_TFR > N_basin);
set(gca,'ydir','normal');
colormap(flipud(gray));
title("Basin ID");

figure;
imagesc((id_Zone_TFR > 0).*(id_Zone_TFR <= N_basin));
set(gca,'ydir','normal');
colormap(flipud(gray));
title("Basin ID");

figure;
imagesc(id_Zone_TFR);
set(gca,'ydir','normal');
colormap(flipud(gray));
title("Basin ID");
pause;

end

