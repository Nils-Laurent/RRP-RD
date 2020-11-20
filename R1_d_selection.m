function [TFR_sig, E_prop_vec] = R1_d_selection(id_TFR, energy_TFR, NB, NZ, P_max)
if exist('P_max','var') == 0
    P_max = 3;
end
E_P_vec = zeros(1, P_max + 1);

[Nfft, L] = size(id_TFR);

% C3_gamma = 3*gamma_Vg;
% isolated_list = id_Zone_TFR.*(id_Zone_TFR <= NB).*(abs(STFT) >= C3_gamma);
% isolated_list = unique(nonzeros(isolated_list(:)));
% E_total = 0;
% for id_zone=1:NB
%     if sum(id_zone == isolated_list) > 0
%         E_total = E_total + Energy_zone(id_zone);
%     end
% end
% for id_zone=(NB + 1):NZ
%     E_total = E_total + Energy_zone(id_zone);
% end
% filt_zones = [isolated_list', (NB + 1):NZ];

filt_TFR = id_TFR > NB;
filt_id_TFR = id_TFR.*filt_TFR;
filt_E_TFR = energy_TFR.*filt_TFR;

E_total = sum(sum(filt_E_TFR));

R_zones = cell(L, 1);
for n=1:L
    [zone_v, zone_k] = sort(energy_TFR(:, n), 'descend');
    if zone_v(1) == 0
        continue;
    end
    [~, s_pos] = unique(nonzeros(zone_v), 'stable');
    zone_k = zone_k(s_pos);
    ids_zones = filt_id_TFR(zone_k, n);
    Rn_zones = nonzeros(ids_zones);
%     r_zones = nonzeros(ids_zones.*(ids_zones > NB));
%     r_zones = intersect(ids_zones, filt_zones, 'stable');
    R_zones{n} = Rn_zones;
end

for P=1:P_max
    
    Z_Pn = zeros(P, L);
    Z_p = zeros(1, L);
    for n=1:L
%         [zone_v, zone_k] = sort(Zone_E_TFR(:, n), 'descend');
%         if zone_v(1) == 0
%             continue;
%         end
%         [~, s_pos] = unique(nonzeros(zone_v), 'stable');
%         zone_k = zone_k(s_pos);
%         ids_zones = id_Zone_TFR(zone_k, n);
%         r_zones = nonzeros(ids_zones.*(ids_zones > NB));
% %         r_zones = intersect(ids_zones, filt_zones, 'stable');

        Rn_zones = R_zones{n, 1};
        if isempty(Rn_zones)
            continue;
        end
        
        P_n = min(P, length(Rn_zones));
        for m=1:P_n
            id_z = Rn_zones(m);
            Z_p(id_z) = id_z;
            Z_Pn(m, n) = id_z;
        end
    end

%     E_P_vec(P + 1) = sum(Energy_zone(Z_p > 0));
    
    TFR_p = zeros(Nfft, L);
    for n=1:L
        for k=1:Nfft
            id_zone = id_TFR(k, n);
            if id_zone > 0 && sum(id_zone == Z_Pn(:, n)) > 0
                TFR_p(k, n) = energy_TFR(k, n);
            end
        end
    end
    E_P_vec(P + 1) = sum(sum(TFR_p));
    

%     figure;
%     imagesc(TFR_p);
%     set(gca,'ydir','normal');
%     colormap(flipud(gray));
%     title("TFR p");
%     pause;
end

TFR_sig = TFR_p;
E_prop_vec = E_P_vec/E_total;

% figure;
% imagesc(TFR_p);
% set(gca,'ydir','normal');
% colormap(flipud(gray));
% title("TFR p");
% 
% figure;
% imagesc(TFR_filt);
% set(gca,'ydir','normal');
% colormap(flipud(gray));
% title("Energy zone filt");

% figure;
% plot((1:(P_max + 1)) - 1, E_P_vec/E_total);
% pause;

return;

%% version avec Delta
% E_vec = zeros(L, 1);
% Delta_vec = zeros(L, 1);
% for n=1:L
%     [zone_v, zone_k] = sort(Zone_E_TFR(:, n), 'ascend');
%     zone_ids = id_Zone_TFR(zone_k, n);
%     id_vec = unique(nonzeros(zone_ids), 'stable');
%     if isempty(id_vec)
%         continue;
%     end
%     E_n = Energy_zone(id_vec)';
%     L_En = length(E_n);
%     E_vec(n) = sum(E_n);
%     Xn = [0, E_n(1:L_En - 1)];
%     delta_n = E_n - Xn;
%     
%     [~, m] = max(delta_n);
%     Delta_vec(n) = L_En - m + 1;
% end
% 
% figure;
% hold on;
% plot(E_vec/max(E_vec));
% plot(Delta_vec, '--');
% hold off;
% title("Energy and delta");
% 
% Delta_energy_vec = zeros(max(Delta_vec), 1);
% for n=1:L
%     m = Delta_vec(n);
%     Delta_energy_vec(m) = Delta_energy_vec(m) + E_vec(n);
% end
% 
% figure;
% plot(Delta_energy_vec);
% title("Delta energy");
% 
% [~, P_est] = max(Delta_energy_vec);
% fprintf("P_est = %u\n", P_est);

% return;
end

