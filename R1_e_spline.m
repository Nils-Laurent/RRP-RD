function [] = R1_e_spline(id_Basins_TFR, id_Zone_TFR, idZ_RP_TFR, E2_Basins_TFR, A_LM_g3, E2_Basins, Energy_Zones, P)

[Nfft, L] = size(id_Zone_TFR);

% delta_s = 3/(sqrt(2*pi)*sigma_s);
% delta_k = ceil(delta_s*Nfft/L);

NB = length(E2_Basins);
NZ = length(Energy_Zones);

filter_TFR = (id_Zone_TFR > NB);
idZF_R_TFR = idZ_RP_TFR.*filter_TFR;
idZF_Z_TFR = id_Zone_TFR.*filter_TFR;
weight_LM_TFR = abs(A_LM_g3.*filter_TFR);
weight_LM_TFR = weight_LM_TFR/max(max(weight_LM_TFR));

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
Modes_max = struct('spline', cell(1, P));
E_max = zeros(1, P);
Mode_Zone_init = zeros(P, NZ);
for n=(n_ord_vec')
    %% set initialization
    if prio_vec(n) == inf
        break;
    end
    
    M_it = struct('spline', cell(1, P));
    acc = 0;
    for m_spl = 1:P
        idZ = KMode_Zone_time(m_spl, n);
        acc = acc + Mode_Zone_init(m_spl, idZ);
        Mode_Zone_init(m_spl, idZ) = 1;
    end
    
    if acc == P
        % initialization did not change
        continue;
    end
    
    %% compute splines
    E_it = zeros(1, P);
    for m_spl = 1:P
        MZ_init = Mode_Zone_init(m_spl, :);
        MZ_other = sum(Mode_Zone_init(1:P ~= m_spl, :), 1);
        %                                   IDZ_R_TFR, ID_Zones_TFR, Energy_TFR, ID_Basins_TFR, Energy_Basins, MZ_init, MZ_other
        [pp, E_it(m_spl)] = R1_e1_spline_it(idZF_R_TFR, idZF_Z_TFR, weight_LM_TFR, id_Basins_TFR, E2_Basins, Energy_Zones, MZ_init, MZ_other);
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
    for m_spl=1:P
        if E_it(m_spl) < E_max(m_spl)
            E_it(m_spl) = E_max(m_spl);
            M_it(m_spl).spline = Modes_max(m_spl).spline;
        end
    end
    
    %% check if modes are crossing
    spl_diff = 0;
    for m_spl=2:P
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
figure;
imagesc((0:L-1)/L, (0:Nfft-1)*L/Nfft, weight_LM_TFR);
set(gca,'ydir','normal');
colormap(flipud(gray));
axis square;
title("spline weight");
hold on;
for m_spl=1:P
    plot((0:L-1)/L, fnval(Modes_max(m_spl).spline, (0:L-1)/L));
end
hold off;

end

