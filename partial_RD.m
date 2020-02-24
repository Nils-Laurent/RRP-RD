function [Weights, k0] = partial_RD(TFR, n0, q, C, gamma, Weights)

[Nfft, L] = size(TFR);

% Values from the Chi-squared table : P[X > Cp] = p
%ratio = 1/100;
%C1 = 9.2103;
ratio = 1/10;
C10 = 4.6052;

TH = gamma*sqrt(C10);

[Magnitue_k0, k0] = max(abs(TFR(:, n0)).*(abs(TFR(:, n0))>3*gamma));
if Magnitue_k0 == 0
    return;
end

Weights(k0, n0) = Weights(k0, n0) + 1;
rq_n0 = round(Nfft/(L^2)*real(q(k0, n0)));

%% forward iterations
rq = rq_n0;
k = k0;
for n=(n0+1):L
    next_index = max(1, min(Nfft, k +rq));
    Im = max(1, next_index -C):min(Nfft, next_index +C);

    Im_ratio = sum(abs(TFR(Im, n)) > TH)/length(Im);
    if Im_ratio <= ratio
        break;
    end

    [~, k] = max(abs(TFR(Im, n)));
    k = k + Im(1)-1;
    Weights(k, n) = Weights(k, n) + 1;
    rq = round(Nfft/(L^2)*real(q(k, n)));
end

%% backward iterations
rq = rq_n0;
k = k0;
for n=(n0-1):-1:1
    next_index = max(1, min(Nfft, k -rq));
    Im = max(1, next_index -C):min(Nfft, next_index +C);

    Im_ratio = sum(abs(TFR(Im, n)) > TH)/length(Im);
    if Im_ratio <= ratio
        break;
    end

    [~, k] = max(abs(TFR(Im, n)));
    k = k + Im(1)-1;
    Weights(k, n) = Weights(k, n) + 1;
    rq = round(Nfft/(L^2)*real(q(k, n)));
end

end

