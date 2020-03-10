function [C_opt] = exridge_n2(TFR, q, sigma_s, gamma_TFR)

[~, L] = size(TFR);
if (nargin < 4)
    gamma_TFR = median(abs(real(TFR(:))))/0.6745;
end

Scale_set = 16:16:128;
N_Scale_Analysis = 3;

NS = length(Scale_set);
curves_1S = zeros(NS, L);
curves_mult_S = zeros(NS-N_Scale_Analysis, L);
energies_mult_S = zeros(NS-N_Scale_Analysis, 1);
E_maxS = 0;

iScale = 0;
for scale=Scale_set
    iScale = iScale + 1;
    fprintf("iScale = %d/%d\n", iScale, NS);
    
    E_max_scale = 0;
    Shift_set=(1 - floor(scale/2)):2:floor(scale/2);
    Base_set = floor(scale/2):scale:(L-(scale/2));
    M_init_n = zeros(length(Shift_set), length(Base_set));
    M_init_k = zeros(length(Shift_set), length(Base_set));
    
    %% scale specific weight computation
    iShift = 0;
        
    Weights = zeros(size(TFR));
    for shift=Shift_set
        iShift = iShift + 1;
        
        sWeights = zeros(size(TFR));

        iN0 = 0;
        for n0=(Base_set + shift)
            iN0 = iN0 + 1;
            [sWeights, k0] = partial_RD(TFR, n0, q, sigma_s, gamma_TFR, sWeights);
            sWeights = sWeights/scale;
            M_init_n(iShift, iN0) = n0;
            M_init_k(iShift, iN0) = k0;
        end
        Weights = max(Weights, sWeights);
    end

    C_opt_scale = zeros(L, 1);
    for n=1:L
        [v, k] = max(sWeights(:, n));
        if v == 0
            continue;
        end

        C_opt_scale(n) = k;
        E_max_scale = E_max_scale + abs(TFR(k, n))^2;
    end
    curves_1S(iScale, :) = C_opt_scale;
    C_opt = C_opt_scale;
    
    %% multi scale analysis
    if iScale > N_Scale_Analysis
        
        for p=N_Scale_Analysis:-1:1
            C_scales = curves_1S(iScale - p, :).*(curves_1S(iScale - p, :) == curves_1S(iScale - p+1, :));
        end

        Energy_S = 0;
        for n=1:L
            if C_scales(n) > 0
                Energy_S = Energy_S + abs(TFR(C_scales(n), n))^2;
            end
        end
        if Energy_S > E_maxS
            E_maxS = Energy_S;
            C_opt = C_scales;
        end
        energies_mult_S(iScale - N_Scale_Analysis) = Energy_S;
        curves_mult_S(iScale - N_Scale_Analysis, :) = C_scales;
    end
end

end