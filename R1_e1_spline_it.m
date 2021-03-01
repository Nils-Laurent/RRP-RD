function [pp_MZ, E_mode] =...
    R1_e1_spline_it(IDZ_R_TFR, ID_Zones_TFR, ID_Basins_TFR,...
    Norm_E_TFR, E_spl_TFR, Energy_Basins, MZ_init, MZ_other, smooth_p, Fs, Nfft)

[N_Y, L] = size(Norm_E_TFR);
% NB = length(Energy_Basins);

TFR_init = zeros(N_Y, L);
for n=1:L
    for k=1:N_Y
        tf_id = ID_Zones_TFR(k, n);
        if tf_id > 0 && MZ_init(tf_id) > 0
            TFR_init(k, n) = Norm_E_TFR(k, n);
        end
    end
end

% figure;
% imagesc((0:L-1)/Fs, (0:N_Y-1)*Fs/Nfft, TFR_init);
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
        k0 = 0;
        w0 = 0;
        for k=1:N_Y
            tf_id = IDZ_R_TFR(k, n);
            if tf_id == 0 || Mode_Zones(tf_id) == 0
                continue;
            end
            if MZ_other(tf_id) > 0
                continue;
            end
            
            Ek = Norm_E_TFR(k, n);
            if Ek > w0
                w0 = Ek;
                k0 = k;
            end
            
        end
        
        if k0 > 0
            m = m + 1;
            iSX(m) = (n - 1)/Fs;
            iSY(m) = (k0 - 1)*Fs/Nfft;
            % iSW(m) = Weight_NB_TFR(k0, n);
            iSW(m) = w0;
        end
    end

    if isempty(iSX)
        break;
    end

    pp_MZ = csaps(iSX, iSY, smooth_p, [], iSW);
    pp_val = fnval(pp_MZ, (0:L-1)/Fs);
    pp_k_vec = round(pp_val*Nfft/Fs) + 1;

%     figure;
%     imagesc((0:L-1)/Fs, (0:N_Y-1)*Fs/Nfft, Energy_TFR);
%     set(gca,'ydir','normal');
%     colormap(flipud(gray));
%     axis square;
%     hold on;
%     plot((0:L-1)/Fs, pp_prev, 'g--');
%     plot((0:L-1)/Fs, pp_val, 'r');
%     plot((0:L-1)/Fs, (pp_k_vec - 1)*Fs/Nfft, 'r:');
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
        if k <= N_Y && k >= 1
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
% fprintf("ITER = %u\n", ITER);
% E_spl_TFR

E_mode = 0;
for n=1:L
    k = pp_k_vec(n);
    if k < 1 || k > N_Y
        continue;
    end
    E_mode = E_mode + E_spl_TFR(k, n);
end

% Mode_Basins = zeros(1, NB);
% for n=1:L
%     k = pp_k_vec(n);
%     if k < 1 || k > N_Y
%         continue;
%     end
%     id_basin = ID_Basins_TFR(k, n);
%     if id_basin > 0
%         Mode_Basins(id_basin) = 1;
%     end
% end
% E_mode = sum(Energy_Basins(Mode_Basins > 0));

% E_mode = sum(Energy_zone(Mode_Zones > 0));

end

