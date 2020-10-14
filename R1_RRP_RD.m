function [] = R1_RRP_RD(STFT, QM, tau, omega, sigma_s, Nr, degree_WPF)

figure;
imagesc(abs(STFT));
set(gca,'ydir','normal');
colormap(flipud(gray));
title("STFT");

%% local maxima TFR
[STFT_LM] = LM_from_STFT(STFT);

g_Vg = median(abs(real(STFT(:))))/0.6745;

Filter_STFT = STFT.*(abs(STFT) > 3*g_Vg);
Filter_LM = STFT_LM.*(abs(STFT) > 3*g_Vg);

%% RRP detection
[id_RRP_TFR, Energy_RRP, ~] = R1_a_idRRP(Filter_LM, STFT, QM);
NB = length(Energy_RRP);

%% Basins of attraction
[id_basin_TFR, Energy_basin] = R1_b_idBasin(STFT, tau, omega, Filter_STFT, STFT_LM, id_RRP_TFR, NB);

[Nfft, L] = size(STFT);

%% Join basins
%DELTA = 5;
% Basins_rel = zeros(NB, 1);
% for n=1:(L-1)
%     for k=1:Nfft
%         basin_id = id_basin_TFR(k, n);
%         if basin_id == 0
%             continue;
%         end
%         
%         k_vec_basin = find(id_basin_TFR(:, n) == basin_id);
%         n_id = id_basin_TFR(k, n);
%         if n_id == 0
%             continue;
%         end
%     end
% end

end

