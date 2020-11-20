function [pp_MZ, E_mode] = R1_e1_spline_it(IDZ_R_TFR, ID_Zones_TFR, Energy_TFR, ID_Basins_TFR, Energy_Basins, Energy_Zones, MZ_init, MZ_other)

[Nfft, L] = size(Energy_TFR);

% TFR_init = zeros(Nfft, L);
% for n=1:L
%     for k=1:Nfft
%         tf_id = ID_TFR(k, n);
%         if tf_id > 0 && MZ_init(tf_id) > 0
%             TFR_init(k, n) = Energy_TFR(k, n);
%         end
%     end
% end
% 
% figure;
% imagesc((0:L-1)/L, (0:Nfft-1)*L/Nfft, TFR_init);
% set(gca,'ydir','normal');
% colormap(flipud(gray));
% axis square;
% title("init");
% pause;
% close all;

pp_prev = zeros(1, L);
MZ_prev = zeros(size(MZ_init));
Mode_Zones = MZ_init;
L_MZ = length(Mode_Zones);
ITER = 0;

while sum(Mode_Zones == MZ_prev) < L_MZ
    ITER = ITER + 1;
    
    m = 0;
    iSX = [];
    iSY = [];
    iSW = [];
    for n=1:L
        for k=1:Nfft
            tf_id = IDZ_R_TFR(k, n);
            if tf_id == 0 || Mode_Zones(tf_id) == 0
                continue;
            end
            if MZ_other(tf_id) > 0
                continue;
            end
            
            m = m + 1;
            iSX(m) = (n - 1)/L;
            iSY(m) = (k - 1)*L/Nfft;
            iSW(m) = Energy_TFR(k, n);
        end
    end
    
    p = 1 - 2/100;
    % spaps (prev. paper)

    pp_MZ = csaps(iSX, iSY, p, [], iSW);
    pp_val = fnval(pp_MZ, (0:L-1)/L);
    pp_k_vec = round(pp_val*Nfft/L) + 1;
    
%     figure;
%     imagesc((0:L-1)/L, (0:Nfft-1)*L/Nfft, Energy_TFR);
%     set(gca,'ydir','normal');
%     colormap(flipud(gray));
%     axis square;
%     hold on;
%     plot((0:L-1)/L, pp_prev, 'g--');
%     plot((0:L-1)/L, pp_val, 'r');
%     hold off;
%     pbaspect([1 1 1]);
%     set(gcf, 'Position',  [0, 0, 1000, 1000]);
%     pause;
%     close all;
%     pp_prev = pp_val;

    MZ_prev = Mode_Zones;
    Mode_Zones = zeros(size(MZ_prev));
    for n=1:L
        k = pp_k_vec(n);
        if k <= Nfft && k >= 1
            tf_id = ID_Zones_TFR(k, n);
            if tf_id > 0
                Mode_Zones(tf_id) = 1;
            end
        end
    end
    
    if ITER == 20
        break;
    end
end
fprintf("ITER = %u\n", ITER);

Mode_Basins = zeros(1, length(Energy_Basins));
for n=1:L
    k = pp_k_vec(n);
    if k < 1 || k > Nfft
        continue;
    end
    id_basin = ID_Basins_TFR(k, n);
    if id_basin > 0
        Mode_Basins(id_basin) = 1;
    end
end
E_mode = sum(Energy_Basins(Mode_Basins > 0));

% E_mode = sum(Energy_zone(Mode_Zones > 0));

end

