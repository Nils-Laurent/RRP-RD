function [id_Basin_TFR, Energy_basin, E2_basin, E_Basins_TFR, E2_Basins_TFR, EB_RP_TFR] =...
    R1_b_idBasin(A_LT_TFR, A_LM_HT_TFR, tau, omega, id_Ridge_TFR, RP_maps, NB, Nfft, Fs)
% A_LT_TFR : |STFT|, low threshold
% A_LM_HT_TFR : Local max of |STFT|, high threshold

[N_Y, L] = size(A_LT_TFR);

id_Basin_TFR = zeros(size(A_LT_TFR));

%% 
NG = 5;
dist_TFR = zeros(size(A_LT_TFR));
for n=(NG+1):(L-NG)
    for k=(NG+1):(N_Y-NG)
        TF_reg = id_Ridge_TFR((k-NG):(k+NG), (n-NG):(n+NG));
        TF_reg = nonzeros(TF_reg);
        if isempty(TF_reg)
            continue;
        end
        dist_TFR(k, n) = mode(TF_reg(:));
    end
end

%% other method : expand LM k to replace zeros in TFR
% need to add STFT_LM as input
% STFT_th = A_LT_TFR
% K_LM_TFR = zeros(size(STFT_th));
% LM_ref = abs(STFT_LM);
% 
% for n=1:L
%     kb = find(LM_ref(:, n), 1);
%     
%     if isempty(kb)
%         continue;
%     else
%         K_LM_TFR(1:kb, n) = kb;
%     end
%     
%     kb_new = kb + find(LM_ref((kb + 1):end, n), 1);
%     while ~(isempty(kb_new))
%         ka = kb;
%         kb = kb_new;
%         k_mid = ceil((kb + ka)/2);
%         
%         K_LM_TFR(ka:k_mid, n) = ka;
%         K_LM_TFR((k_mid+1):kb, n) = kb;
%         
%         kb_new = kb + find(LM_ref((kb + 1):end, n), 1);
%     end
%     
%     K_LM_TFR(kb:end, n) = kb;
% end

% WHEN IN THE LOOP: 
%         w_LM = K_LM_TFR(w, t);
%         if w_LM == 0
%             continue;
%         end
%         
%         basin_id = id_Ridge_TFR(w_LM, t);

%% save vectors for all coefficients
RT_g = round(real(tau)*Fs);
for n=1:L
    RT_g(:, n) = RT_g(:, n) + n;
end

RT_g(RT_g > L) = L;
RT_g(RT_g < 1) = 1;

Omega_g = round(omega*Nfft/Fs) + 1;

Omega_g(Omega_g > N_Y) = N_Y;
Omega_g(Omega_g < 1) = 1;

Energy_basin = zeros(NB, 1);
%% Create basins allowing detached coefficients
for n=1:L
    for k=1:N_Y
        if A_LT_TFR(k, n) == 0
            continue;
        end
        
        t = RT_g(k, n);
        w = Omega_g(k, n);
        
        basin_id = dist_TFR(w, t);
        if basin_id == 0
            continue;
        end
        
        %% add coefficient to basin
        Energy_basin(basin_id) = Energy_basin(basin_id) + abs(A_LT_TFR(k, n));
        id_Basin_TFR(k, n) = basin_id;
    end
end


% figure;
% imagesc(id_Basin_TFR);
% set(gca,'ydir','normal');
% colormap(flipud(gray));
% title("id basins");

%% remove basins without important ridges
E2_basin = zeros(1, NB);
map_rm = ones(1, NB);

for n=1:L
    for k=1:N_Y
        e_TF = A_LM_HT_TFR(k, n);
        if e_TF == 0
            continue;
        end
        
        id_ridge = id_Ridge_TFR(k, n);
        if id_ridge > 0 && id_ridge == id_Basin_TFR(k, n)
            E2_basin(id_ridge) = E2_basin(id_ridge) + e_TF;
            map_rm(id_ridge) = 0;
        end
    end
end

% for n=1:L
%     for k=1:Nfft
%         id_ridge = id_Ridge_TFR(k, n);
%         if id_ridge == 0
%             continue;
%         end
%         
%         if id_ridge == id_Basin_TFR(k, n)
%             map_rm(id_ridge) = 0;
%         end
%     end
% end

for n=1:L
    for k=1:N_Y
        id_basin = id_Basin_TFR(k, n);
        if id_basin == 0
            continue;
        end
        if map_rm(id_basin) == 1
            Energy_basin(id_basin) = 0;
            id_Basin_TFR(k, n) = 0;
        end
    end
end


%% remove detached coefficients
id_Basin_TFR2 = id_Basin_TFR;
for n0=1:L
    for k0=1:N_Y
        basin_id = id_Basin_TFR(k0, n0);
        if basin_id == 0
            continue;
        end
        
        n = RT_g(k0, n0);
        k = Omega_g(k0, n0);
        
        dn = n - n0;
        dk = k - k0;
        sk = sign(dk) + (dk == 0);
        
        if dn == 0 && sum(abs(id_Basin_TFR(k0:sk:k, n) - basin_id)) > 0
            %% remove coefficient from basin
            id_Basin_TFR2(k0, n0) = 0;
            Energy_basin(basin_id) = Energy_basin(basin_id) - abs(A_LT_TFR(k0, n0));
            continue;
        end
        
        %% compute line
        [n_vec, k_vec] = bresenham(n0, k0, n, k);
        max_tf = 0;
        for m=1:length(n_vec)
            ni = n_vec(m);
            ki = k_vec(m);
            max_tf = max_tf + abs(id_Basin_TFR(ki, ni) - basin_id);
            if max_tf > 0
                break;
            end
        end
        
        if max_tf > 0
            %% remove coefficient from basin
            id_Basin_TFR2(k0, n0) = 0;
            Energy_basin(basin_id) = Energy_basin(basin_id) - abs(A_LT_TFR(k0, n0));
            continue;
        end
    end
end
id_Basin_TFR = id_Basin_TFR2;

E_Basins_TFR = zeros(size(A_LT_TFR));
E2_Basins_TFR = zeros(size(A_LT_TFR));
%% basin energy
for n=1:L
    for k=1:N_Y
        basin_id = id_Basin_TFR(k, n);
        if basin_id > 0
            E_Basins_TFR(k, n) = Energy_basin(basin_id);
            E2_Basins_TFR(k, n) = E2_basin(basin_id);
        end
    end
end

EB_RP_TFR = zeros(size(A_LT_TFR));
for m=1:RP_maps.M_map
    RE = Energy_basin(RP_maps.ids(m));
    n = RP_maps.ns(m);
    k = RP_maps.ks(m);
    EB_RP_TFR(k, n) = RE;
end

end

