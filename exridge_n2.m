function [C_opt, energies] = exridge_n2(TFR, q, C, R_set)

[Nfft, L] = size(TFR);

gamma = median(abs(real(TFR(:))))/0.6745;

% Values from the Chi-squared table : P[X > Cp] = p
%ratio = 1/100;
%C1 = 9.2103;
ratio = 1/10;
C10 = 4.6052;

TH = gamma*sqrt(C10);

E_max = 0;
%R_set = 8:8:512;
NR = length(R_set);
energies = zeros(NR, 1);
curves = zeros(NR, L);
ind_R = 0;

for shift_size=R_set
    ind_R = ind_R + 1;
    if mod(ind_R, 16) == 0
        fprintf("ind_R = %d/%d\n", ind_R, NR);
    end
    Shifts = 1:shift_size:L;
    N_shift = length(Shifts);

    Weights = zeros(size(TFR));
    Step_mat = zeros(N_shift, L);

    for n0=Shifts
        [e0, y0] = max(abs(TFR(:, n0)).*(abs(TFR(:, n0))>3*gamma));
        if e0 == 0
            continue;
        end

        Weights(y0, n0) = Weights(y0, n0) + 1;
        Step_mat(y0, n0) = n0;
        rq_n0 = round(Nfft/(L^2)*real(q(y0, n0)));

        %% forward iterations
        rq = rq_n0;
        arg = y0;
        for n=(n0+1):L
            next_index = min(Nfft, arg +rq);
            next_index = max(1, next_index);
            Im = max(1, next_index -C):min(Nfft, next_index +C);

            Im_ratio = sum(abs(TFR(Im, n)) > TH)/length(Im);
            if Im_ratio <= ratio
                break;
            end

            [~, arg] = max(abs(TFR(Im, n)));
            arg = arg + Im(1)-1;
            Weights(arg, n) = Weights(arg, n) + 1;
            Step_mat(arg, n) = n0;
            rq = round(Nfft/(L^2)*real(q(arg, n)));
        end

        %% backward iterations
        rq = rq_n0;
        arg = y0;
        for n=(n0-1):-1:1
            next_index = min(Nfft, arg -rq);
            next_index = max(1, next_index);
            Im = max(1, next_index -C):min(Nfft, next_index +C);

            Im_ratio = sum(abs(TFR(Im, n)) > TH)/length(Im);
            if Im_ratio <= ratio
                break;
            end

            [~, arg] = max(abs(TFR(Im, n)));
            arg = arg + Im(1)-1;
            Weights(arg, n) = Weights(arg, n) + 1;
            Step_mat(arg, n) = n0;
            rq = round(Nfft/(L^2)*real(q(arg, n)));
        end
    end

    C_test = zeros(L, 1);
    E_test = 0;
    for n=1:L
        [v, k] = max(Weights(:, n));
        if v == 0
            continue;
        end
        C_test(n) = k;
        curves(ind_R, n) = k;
        E_test = E_test + abs(TFR(k, n))^2;
        
    end
    energies(ind_R) = E_test;
    if E_test > E_max
        E_max = E_test;
        C_opt = C_test;
    end
end

% figure;
% plot(energies);

end