function [id_Basin_TFR, Energy_basin, Energy_TFR, E_RB_TFR] =...
    R1_b_idBasin(STFT_th, tau, omega, STFT_LM, id_Ridge_TFR, RP_maps, NB)

[Nfft, L] = size(STFT_th);

id_Basin_TFR = zeros(size(STFT_th));
Energy_basin = zeros(NB, 1);
Energy_TFR = zeros(size(STFT_th));
E_RB_TFR = zeros(size(STFT_th));

%% test
K_LM_TFR = zeros(size(STFT_th));
LM_ref = abs(STFT_LM);

for n=1:L
    kb = find(LM_ref(:, n), 1);
    
    if isempty(kb)
        continue;
    else
        K_LM_TFR(1:kb, n) = kb;
    end
    
    kb_new = kb + find(LM_ref((kb + 1):end, n), 1);
    while ~(isempty(kb_new))
        ka = kb;
        kb = kb_new;
        k_mid = ceil((kb + ka)/2);
        
        K_LM_TFR(ka:k_mid, n) = ka;
        K_LM_TFR((k_mid+1):kb, n) = kb;
        
        kb_new = kb + find(LM_ref((kb + 1):end, n), 1);
    end
    
    K_LM_TFR(kb:end, n) = kb;
end

% figure;
% imagesc(K_LM_TFR);
% set(gca,'ydir','normal');
% colormap(flipud(gray));
% title("K LM");
% pause;

RT_g = round(real(tau)*L);
for n=1:L
    RT_g(:, n) = RT_g(:, n) + n;
end

RT_g(RT_g > L) = L;
RT_g(RT_g < 1) = 1;

Omega_g = round(omega*Nfft/L) + 1;

Omega_g(Omega_g > Nfft) = Nfft;
Omega_g(Omega_g < 1) = 1;

%% Basins TFR code
for n=1:L
    for k=1:Nfft
        if STFT_th(k, n) == 0
            continue;
        end
        
        t = RT_g(k, n);
        w = Omega_g(k, n);
        
        w_LM = K_LM_TFR(w, t);
        if w_LM == 0
            continue;
        end
        
        basin_id = id_Ridge_TFR(w_LM, t);
        if basin_id == 0
            continue;
        end
        Energy_basin(basin_id) = Energy_basin(basin_id) + abs(STFT_th(k, n))^2;
        id_Basin_TFR(k, n) = basin_id;
    end
end

% RP_maps.ids = id_map_RRP;
% RP_maps.ns = n_map_RRP;
% RP_maps.ks = k_map_RRP;
% RP_maps.M_map = M_map;
for m=1:RP_maps.M_map
    RE = Energy_basin(RP_maps.ids(m));
    n = RP_maps.ns(m);
    k = RP_maps.ks(m);
    E_RB_TFR(k, n) = RE;
end

% for basin_id=1:NB
% %     Energy_TFR = Energy_TFR + Energy_basin(basin_id)*(id_Basin_TFR == basin_id);
%     E_RB_TFR = E_RB_TFR + Energy_basin(basin_id)*(id_Ridge_TFR == basin_id);
% end


% figure;
% imagesc(Energy_TFR);
% set(gca,'ydir','normal');
% colormap(flipud(gray));
% title("Basins energy");
end

