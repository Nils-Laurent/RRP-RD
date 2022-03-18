function [R_out, E_out] = RRP_RD(STFT, Sc, omega, tau, P, varargin)
% RRP_RD extracts the ridges of a multicomponent signal
% [R_out, E_out] = RRP_RD(STFT, Sc, omega, tau, P)
% [R_out, E_out] = RRP_RD(STFT, Sc, omega, tau, P, Nr)
%
%   Extracts the ridge curve by maximising some energy.
%   The algorithm uses relevant ridge portions to detect the ridge [1].
%
% INPUTS:
%   STFT    : short time Fourier transform.
%   Sc      : Chirp Rate (CR) estimate [2].
%   omega   : First order Instantaneous Frequency (IF) estimate [2].
%   tau     : groupe delay [2].
%   P       : spline smoothing parameter,
%             typicall value for a high level of noise is 1 - 1E-4,
%                                  low level of noise is 1 - 1E-6.
%
%   Nr      : number of ridges
%             equals to 1 by default.
%
%   --  List of possible name-value pair argument
%   'samp'  : samples per second,
%             equals to size(STFT, 2) by default.
%   'nfft'  : number of frequency bins,
%             equals to size(STFT, 1) by default.
%
% OUTPUTS:
%   R_out   : array of length Nr containing structures with field 'spline',
%             for example, to get the spline corresponding to the n-th
%             ridge, the code would be 'R_out(n).spline'.
%
%   E_out   : array of length Nr containing the energy associated with the
%             detection of each ridge in R_out.
%
% REFERENCES:
% [1] N. Laurent and S. Meignen, "A Novel Ridge Detector for Nonstationary
%     Multicomponent Signals: Development and Application to Robust Mode
%     Retrieval," in IEEE Transactions on Signal Processing, vol. 69,
%     pp. 3325-3336, 2021, doi: 10.1109/TSP.2021.3085113.
%
% [2] Behera, R., Meignen, S., & Oberlin, T. (2015). Theoretical Analysis
%     of the Second-order Synchrosqueezing Transform. To appear in ACHA

defaultSamp = size(STFT, 2);
defaultNfft = size(STFT, 1);
defaultNr = 1;

p = inputParser;
addRequired(p,'STFT');
addRequired(p,'Sc');
addRequired(p,'omega');
addRequired(p,'tau');
addRequired(p,'P');
addOptional(p,'Nr',defaultNr);
addParameter(p,'samp',defaultSamp);
addParameter(p,'nfft',defaultNfft);
parse(p,STFT,Sc,omega,tau,P,varargin{:});

Nr = p.Results.Nr;
samp = p.Results.samp;
Nfft = p.Results.nfft;

[STFT_LM] = LM_from_STFT(STFT);

gamma_Vg = median(abs(real(STFT(:))))/0.6745;

C2_gamma = 2*gamma_Vg;
C3_gamma = 3*gamma_Vg;

ASTFT_g2 = abs(STFT).*(abs(STFT) > C2_gamma);
A_LM_g2 = abs(STFT_LM).*(abs(STFT) > C2_gamma);
A_LM_g3 = abs(STFT_LM).*(abs(STFT) > C3_gamma);

[id_RP_TFR, Energy_RP, E_RP_TFR, RP_maps] =...
    R1_a_idRRP(A_LM_g2, STFT, Sc, Nfft, samp);
NB = length(Energy_RP);

% figure;
% imagesc(A_LM_g2);
% set(gca,'ydir','normal');
% colormap(flipud(gray));
% title("LM g2");
% 
% figure;
% imagesc(A_LM_g3);
% set(gca,'ydir','normal');
% colormap(flipud(gray));
% title("LM g3");
% pause;

[id_Basins_TFR, Energy_basins, E2_basins, E_Basins_TFR, E2_Basins_TFR, EB_RP_TFR] =...
    R1_b_idBasin(ASTFT_g2, A_LM_g3, tau, omega, id_RP_TFR, RP_maps, NB, Nfft, samp);

% figure;
% imagesc(E_Basins_TFR);
% set(gca,'ydir','normal');
% colormap(flipud(gray));
% title("E basins");

[id_Zones_TFR, Energy_Zones, idZ_RP_TFR, EZ_RP_TFR, E_Zones_TFR] =...
    R1_c_idZones(ASTFT_g2, id_Basins_TFR, Energy_basins, id_RP_TFR);

% figure;
% imagesc(E_Zones_TFR);
% set(gca,'ydir','normal');
% colormap(flipud(gray));
% title("E zones");
% pause;

% fprintf("Set splines\n");

[R_out, E_out] = R1_e_spline(E_Zones_TFR, id_Basins_TFR, id_Zones_TFR,...
    idZ_RP_TFR, E2_Basins_TFR, A_LM_g3, E2_basins, Energy_Zones, Nr, P, samp, Nfft);

return;

%% analyse des zones

% g_basin_disp = id_Basin_TFR + 14*(id_Basin_TFR > 0);
% figure;
% imagesc(g_basin_disp);
% set(gca,'ydir','normal');
% % myColorMap = colorcube(round(NB/8));
% myColorMap = colorcube;
% myColorMap(1,:) = 1;
% colormap(myColorMap);
% title("Basin ID");
% 
% g_zone_disp = id_Zone_TFR.*(id_Zone_TFR > NB);
% g_zone_disp = g_zone_disp - NB*(id_Zone_TFR > NB);
% 
% figure;
% imagesc(g_zone_disp);
% set(gca,'ydir','normal');
% myColorMap = lines(NZ - NB);
% myColorMap(1,:) = 1;
% colormap(myColorMap);
% title("Zone ID");
% return;

end

