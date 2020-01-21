function [c] = exridge_variant(TFR, sigma, q, C)


[Nfft, L] = size(TFR);
c = zeros(L, 1);


TFR_it = TFR;

[~, arg] = max(abs(TFR_it(:)));
[y0, n0] = ind2sub(size(TFR),arg);
c(n0) = y0;

RQ = zeros(L, 1);
RQ(n0) = round(Nfft/(L^2)*real(q(c(n0), n0)))+1;

%% forward iterations
for m=(n0+1):L
    Im = max(1, c(m-1) +RQ(m-1) -C):min(Nfft, c(m-1) +RQ(m-1) +C);
    [~, arg] = max(abs(TFR_it(Im, m)));
    c(m) = arg + Im(1)-1;
    RQ(m) = round(Nfft/(L^2)*real(q(c(m), m)))+1;
end

%% backward iterations
for m=(n0-1):-1:1
    Im = max(1, c(m+1) +RQ(m+1) -C):min(Nfft, c(m+1) +RQ(m+1) +C);
    [~, arg] = max(abs(TFR_it(Im, m)));
    c(m) = arg + Im(1)-1;
    RQ(m) = round(Nfft/(L^2)*real(q(c(m), m)))+1;
end

for n=1:L
    B = real(q(c(n), n));
    eta_k = ceil(Nfft/L*3/sqrt(2*pi)*sqrt(1/sigma^2 + sigma^2*B^2));
    range = max(1, c(n)-eta_k):min(Nfft, c(n)+eta_k);
    TFR_it(range, n) = 0;
end

figure;
imagesc((0:L-1)/L, (L/Nfft)*(1:Nfft), abs(TFR_it));
set(gca,'ydir','normal');
title("debug");
axis square
hold on;
plot((0:L-1)/L, c*L/Nfft, 'r');
hold off;

end