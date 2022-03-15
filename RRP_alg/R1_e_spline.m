function [Modes_max, E_max] = R1_e_spline(E_Zones_TFR, id_Basins_TFR, id_Zones_TFR,...
    idZ_RP_TFR, E2_Basins_TFR, A_LM_g3, E2_Basins, Energy_Zones, Nr, smooth_p, Fs, Nfft)

[N_Y, L] = size(id_Zones_TFR);

% delta_s = 3/(sqrt(2*pi)*sigma_s);
% delta_k = ceil(delta_s*Nfft/L);

NB = length(E2_Basins);
NZ = length(Energy_Zones);

%% prepare TFR for mode energy computation
E_spl_TFR = zeros(N_Y, L);
for n=1:L
    EB_n = zeros(1, NB);
    for k=1:N_Y
        E_kn = A_LM_g3(k, n);
        if E_kn == 0
            continue;
        end
        
        id_basin = id_Basins_TFR(k, n);
        if id_basin > 0
            EB_n(id_basin) = E_kn;
        end
    end
    
    for k=1:N_Y
        id_basin = id_Basins_TFR(k, n);
        if id_basin > 0
            E_spl_TFR(k, n) = EB_n(id_basin);
        end
    end
end

% filter_TFR = (id_Zones_TFR > 0);
% idZF_R_TFR = idZ_RP_TFR.*filter_TFR;
% idZF_Z_TFR = id_Zones_TFR.*filter_TFR;
% EZF_TFR = E_Zones_TFR.*filter_TFR;

%% create weight matrices
Norm_g3 = A_LM_g3;
% Norm_g3 = Norm_g3/max(max(Norm_g3));

% Weight_NB_TFR = Norm_g3;
% Div_Basins = zeros(1, NB);
% for n=1:L
%     for k=1:N_Y
%         id_basin = id_Basins_TFR(k, n);
%         if id_basin > 0
%             if Norm_g3(k, n) > Div_Basins(id_basin)
%                 Div_Basins(id_basin) = Norm_g3(k, n);
%             end
%         end
%     end
% end
% for n=1:L
%     for k=1:N_Y
%         id_basin = id_Basins_TFR(k, n);
%         if id_basin > 0 && Div_Basins(id_basin) > 0
%             Weight_NB_TFR(k, n) = Weight_NB_TFR(k, n)/Div_Basins(id_basin);
%         end
%     end
% end

%% used TFR
U_E_TFR = E_Zones_TFR;
U_ID_TFR = id_Zones_TFR;
U_IDR_TFR = idZ_RP_TFR;
% idZF_R_TFR = idZ_RP_TFR.*filter_TFR;
% idZF_Z_TFR = id_Zones_TFR.*filter_TFR;
% EZF_TFR = E_Zones_TFR.*filter_TFR;

%% new init

ID_vec = zeros(L, Nr);
for n=1:L
    E_vec = U_E_TFR(:, n);
    [ES_vec, sort_k] = sort(E_vec, 'descend');
    [~, u_pos] = unique(nonzeros(ES_vec), 'stable');
    if length(u_pos) >= Nr
        k_vec = sort_k(u_pos(1:Nr));
        k_vec = sort(k_vec, 'ascend');
        id_vec = U_ID_TFR(k_vec, n);
        if min(id_vec) == 0
            continue;
        end
        ID_vec(n, :) = id_vec;
    end
end

ID_vec = sortrows(ID_vec);
occ_vec = zeros(1, L);
if sum(ID_vec(1, :)) > 0
    occ_vec(1) = 1;
end
for n=2:L
    if sum(ID_vec(n, :)) == 0
        continue;
    end
    if sum(ID_vec(n, :) == ID_vec(n - 1, :)) == Nr
        occ_vec(n) = occ_vec(n - 1) + 1;
        occ_vec(n - 1) = 0;
    else
        occ_vec(n) = 1;
    end
end

[~, n_ord_vec1] = sort(occ_vec, 'descend');
KMode_Zone_time = ID_vec.';
n_ord_vec = n_ord_vec1.';

%% splines
Modes_max = struct('spline', cell(1, Nr));
E_max = zeros(1, Nr);
Mode_Zone_init = zeros(Nr, NZ);
for n=(n_ord_vec')
    %% set initialization
%     if prio_vec(n) == inf
%         break;
%     end
    if occ_vec(n) == 0
        break;
    end
    
    M_it = struct('spline', cell(1, Nr));
    acc = 0;
    for m_spl = 1:Nr
        idZ = KMode_Zone_time(m_spl, n);
        acc = acc + Mode_Zone_init(m_spl, idZ);
        Mode_Zone_init(m_spl, idZ) = 1;
    end
    
    if acc == Nr
        % initialization did not change
        continue;
    end
    
    %% compute splines
    E_it = zeros(1, Nr);
    for m_spl = 1:Nr
        MZ_init = Mode_Zone_init(m_spl, :);
        MZ_other = sum(Mode_Zone_init(1:Nr ~= m_spl, :), 1);
        [pp, E_it(m_spl)] = R1_e1_spline_it(U_IDR_TFR, U_ID_TFR, id_Basins_TFR,...
            Norm_g3, E_spl_TFR, E2_Basins, MZ_init, MZ_other, smooth_p, Fs, Nfft);
        M_it(m_spl).spline = pp;
    end
    
    %% figure
%     E_it
%     figure;
%     imagesc((0:L-1)/L, (0:Nfft-1)*L/Nfft, E2_Basins_TFR);
%     set(gca,'ydir','normal');
%     colormap(flipud(gray));
%     axis square;
%     title("spline weight");
%     hold on;
%     for m_spl=1:P
%         plot((0:L-1)/L, fnval(M_it(m_spl).spline, (0:L-1)/L));
%     end
%     hold off;
%     pause;
%     close all;
    
    %% keep modes with maximum energy
    for m_spl=1:Nr
        if E_it(m_spl) < E_max(m_spl)
            E_it(m_spl) = E_max(m_spl);
            M_it(m_spl).spline = Modes_max(m_spl).spline;
        end
    end
    
    %% check if modes are crossing
    spl_diff = 0;
    for m_spl=2:Nr
        y_hf = fnval(M_it(m_spl).spline, (0:L-1)/L);
        y_lf = fnval(M_it(m_spl - 1).spline, (0:L-1)/L);
        spl_diff = min(y_hf - y_lf);
        if spl_diff < 0
            break;
        end
    end
    if spl_diff < 0
        % modes are crossing
        continue;
    end
    
    E_max = E_it;
    Modes_max = M_it;
end

    
%% figure
% figure;
% imagesc((0:L-1)/L, (0:N_Y-1)*Fs/Nfft, Weight_NB_TFR);
% set(gca,'ydir','normal');
% colormap(flipud(gray));
% axis square;
% title("spline weight");
% hold on;
% for m_spl=1:Nr
%     plot((0:L-1)/Fs, fnval(Modes_max(m_spl).spline, (0:L-1)/Fs));
% end
% hold off;
% pause;
% close all;

% X = fnval(fnder(fnder(Modes_max(1).spline)), (0:L-1)/L);
% Y = 1/L*norm(X)^2

return;

%% old init
[~, idZ_vec] = sort(Energy_Zones, 'descend');
[~, spZ_vec] = sort(idZ_vec);

KMode_Zone_time = zeros(P, L);
prio_vec = inf(L, P);
for n=1:L
%     idVec = idZF_R_TFR(:, n);
%     spVec = zeros(size(idVec));
%     %% set delta separation vec
%     for k=1:Nfft
%         id_zone = idVec(k);
%         if id_zone > 0
%             spVec(k) = spZ_vec(id_zone);
%         end
%     end
%     
%     for k=1:Nfft
%         if spVec(k) == 0
%             continue;
%         end
%         
%         ka = max(1, k - delta_k);
%         kb = min(Nfft, k + delta_k);
%         min_v = min(nonzeros(spVec(ka:kb)));
%         spVec(ka:kb) = min_v*(spVec(ka:kb) == min_v);
%     end
%     
%     idVec = idVec.*(spVec > 0);

%     idVec = idZF_R_TFR(:, n);
%     %% EB new TFR
    idVec = idZF_R_TFR(:, n).*(E2_Basins_TFR(:, n) > 0);
    
    %% find modes at n
    SPn = inf(1, Nfft);
    for k=1:Nfft
        % idZ = idF_R_TFR(k, n);
        idZ = idVec(k);
        if idZ == 0
            continue;
        end
        
        SPn(k) = spZ_vec(idZ);
    end
    [SPn, s_pos] = sort(SPn, 'ascend');
    [SPn, u_pos] = unique(SPn(SPn < inf), 'stable');
    
    %% set priority vec and frequencies
    if length(SPn) >= P
        %% P_XX are the mode's freq. bin in freq. order
        [~, P_XX] = sort(s_pos(u_pos(1:P)));
        for m=1:P
            prio_vec(n, m) = SPn(P - m + 1);
            KMode_Zone_time(P_XX(m), n) = idZ_vec(SPn(m));
        end
    end
end

% [~, n_ord_vec] = sort(prio_vec, 'ascend');
[~, n_ord_vec] = sortrows(prio_vec);
end

