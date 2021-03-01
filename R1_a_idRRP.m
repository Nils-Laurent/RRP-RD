function [id_TFR, Energy_RRP, Energy_TFR, RP_maps] = R1_a_idRRP(STFT_LM_th, STFT, QM, Nfft, Fs)

ASTFT_LM_th = abs(STFT_LM_th);

% figure;
% imagesc(ASTFT_LM_th > 0);
% set(gca,'ydir','normal');
% colormap(flipud(gray));
% title("LM th 1");
% pause;

N_BS_RRP = sum(sum(ASTFT_LM_th > 0));
Energy_RRP = zeros(N_BS_RRP, 1);

id_map_RRP = zeros(N_BS_RRP, 1);
k_map_RRP = zeros(N_BS_RRP, 1);
n_map_RRP = zeros(N_BS_RRP, 1);
M_map = 0;

Energy_TFR = zeros(size(STFT));

%% identify RRP

[N_Y, L] = size(ASTFT_LM_th);

N_RRP = 0;
id_TFR = zeros(size(ASTFT_LM_th));
rSPEC = abs(STFT);
for n=1:(L-1)
    [vec_LM] = find(ASTFT_LM_th(:, n));
    [next_vec_LM] = find(ASTFT_LM_th(:, n+1));
    
    if isempty(next_vec_LM) || isempty(vec_LM)
        continue;
    end
    
    for k=2:(N_Y-1)
        if ASTFT_LM_th(k, n) == 0
            continue;
        end
        
        %% next
        rqa_grid = round(real(QM(k, n))*Nfft/(Fs^2));
        kb = max(1, min(N_Y, k + rqa_grid));
        
        [~, arg_w] = min(abs(next_vec_LM - kb));
        kb_LM = next_vec_LM(arg_w);
        
        %% next_prev
        rqb_grid = round(real(QM(kb_LM, n+1))*Nfft/(Fs^2));
        ka = max(1, min(N_Y, kb_LM - rqb_grid));
        
        [~, arg_ka] = min(abs(vec_LM - ka));
        ka_LM = vec_LM(arg_ka);
        
        %% equivalent
        if (k == ka_LM)
            if id_TFR(k, n) > 0
                id_RRP = id_TFR(k, n);
            else
                % New RRP
                N_RRP = N_RRP + 1;
                id_RRP = N_RRP;
                id_TFR(k, n) = id_RRP;
                Energy_RRP(id_RRP) = Energy_RRP(id_RRP) + rSPEC(k, n);
                
                % map RRP
                M_map = M_map + 1;
                id_map_RRP(M_map) = id_RRP;
                n_map_RRP(M_map) = n;
                k_map_RRP(M_map) = k;
            end
            id_TFR(kb_LM, n+1) = id_RRP;
            Energy_RRP(id_RRP) = Energy_RRP(id_RRP) + rSPEC(kb_LM, n+1);
            
            % map RRP
            M_map = M_map + 1;
            id_map_RRP(M_map) = id_RRP;
            n_map_RRP(M_map) = n+1;
            k_map_RRP(M_map) = kb_LM;
        end
    end
end
Energy_RRP = Energy_RRP(1:N_RRP);

id_map_RRP = id_map_RRP(1:M_map);
n_map_RRP = n_map_RRP(1:M_map);
k_map_RRP = k_map_RRP(1:M_map);

%% compute RRP energy
for m=1:M_map
    RE = Energy_RRP(id_map_RRP(m));
    n = n_map_RRP(m);
    k = k_map_RRP(m);
    Energy_TFR(k, n) = RE;
end

RP_maps.ids = id_map_RRP;
RP_maps.ns = n_map_RRP;
RP_maps.ks = k_map_RRP;
RP_maps.M_map = M_map;
% figure;
% imagesc(Energy_TFR);
% set(gca,'ydir','normal');
% colormap(flipud(gray));
% title("Energy TFR");
% pause;

end

