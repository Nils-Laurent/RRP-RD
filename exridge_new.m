function [Cr] = exridge_new(TFR, sigma, q, C)

[Nfft, L] = size(TFR);
%c = zeros(L, 1);

TFR_it = TFR;

gamma = median(abs(real(TFR(:))))/0.6745;

% Values from the Chi-squared table : P[X > Cp] = p
%C1 = 9.2103;
ratio = 1/10;
C10 = 4.6052;

TH = gamma*sqrt(C10);

Nr = min(60,floor(L/8));
init_set = floor(linspace(L/((Nr+1)), L-L/((Nr+1)), Nr));
init_off = 0:(Nr-1);
r = 0;

Cr = zeros(Nr, Nr, L);
Ar = zeros(Nr, 1);
Br = zeros(Nr, 1);
Er = zeros(Nr, 1);

X_init = zeros(Nr, 1);
Y_init = zeros(Nr, 1);

%for c_off = init_off
c_off = 0;
for nd = (init_set + c_off)
    r = r + 1;
    [e, y0] = max(abs(TFR_it(:, nd)));
    X_init(r) = nd;
    Y_init(r) = y0;
    e = e^2;
    Cr(r, nd) = y0;
    Ar(r) = nd;
    Br(r) = nd;

    RQ = zeros(L, 1);
    %RQ = RQ + round(4000*(Nfft/(L^2)));
    RQ(nd) = round(Nfft/(L^2)*real(q(Cr(r, nd), nd)));

    %% forward iterations
    for N_parts=(nd+1):L
        Im = max(1, Cr(r, N_parts-1) +RQ(N_parts-1) -C):min(Nfft, Cr(r, N_parts-1) +RQ(N_parts-1) +C);
        
        Im_ratio = sum(abs(TFR_it(Im, N_parts)) > TH)/length(Im);
        if Im_ratio <= ratio
            break;
        end

        [v_arg, arg] = max(abs(TFR_it(Im, N_parts)));
        e = e + v_arg^2;
        Br(r) = N_parts;
        Cr(r, N_parts) = arg + Im(1)-1;
        RQ(N_parts) = round(Nfft/(L^2)*real(q(Cr(r, N_parts), N_parts)));
    end

    %% backward iterations
    for N_parts=(nd-1):-1:1
        Im = max(1, Cr(r, N_parts+1) -RQ(N_parts+1) -C):min(Nfft, Cr(r, N_parts+1) -RQ(N_parts+1) +C);

        Im_ratio = sum(abs(TFR_it(Im, N_parts)) > TH)/length(Im);
        if Im_ratio <= ratio
            break;
        end

        [v_arg, arg] = max(abs(TFR_it(Im, N_parts)));
        e = e + v_arg^2;
        Ar(r) = N_parts;
        Cr(r, N_parts) = arg + Im(1)-1;
        RQ(N_parts) = round(Nfft/(L^2)*real(q(Cr(r, N_parts), N_parts)));
    end
    
    Er(r) = e;
end
%end

[Er, Er_order] = sort(Er, 'descend');
Br = Br(Er_order);
Ar = Ar(Er_order);
Cr = Cr(Er_order, :);
X_init = X_init(Er_order);
Y_init = Y_init(Er_order);

%% Ridge assignement
opt_map = zeros(L, 1);
C_opt = zeros(L, 1);
Nr_opt = 0;
Er_opt = zeros(Nr, 2);
for r=1:Nr
    nd = X_init(r);
    if opt_map(nd) > 0
        continue;
    end
    
    Nr_opt = Nr_opt + 1;
    Er_opt(Nr_opt, 2) = r;
    len_opt_r = 0;
    for m = Ar(r):Br(r)
        if opt_map(m) == 0
            len_opt_r = len_opt_r + 1;
            opt_map(m) = r;
            C_opt(m) = Cr(opt_map(m), m);
            Er_opt(Nr_opt, 1) = Er_opt(Nr_opt, 1) + abs(TFR(C_opt(m), m))^2;
        end
    end
    
    if Er_opt(Nr_opt, 1)/len_opt_r < TH^2
        fprintf("low energy ridge = %d\n", Er_opt(Nr_opt, 1));
    end
end
C_opt(C_opt == 0) = nan;

QR_opt = zeros(L, 1);
for m=init_set
    if ~isnan(C_opt(m))
        QR_opt(m) = real(q(C_opt(m), m));
    else
        QR_opt(m) = nan;
    end
end

% figure;
% plot(QR_opt);
% pause;

% figure;
% plot(Er_opt(:, 1));
% pause;

%% keep ridges where initialization is used
% for n=1:L
%     if opt_map(n) > 0
%         n0 = X_init1(opt_map(n));
%         if opt_map(n0) ~= opt_map(n)
%             opt_map(opt_map == opt_map(n)) = 0;
%         end
%     end
% end

N_parts = 0;
E_parts = [];
X_delim = [];
Y_delim = [];
X_init_opt = [];
Y_init_opt = [];
for m=2:L
    if opt_map(m-1) == opt_map(m)
        continue;
    end
    if opt_map(m) == 0
        continue;
    end
    
    N_parts = N_parts + 1;
    r = opt_map(m);
    E_parts(N_parts) = Er_opt(Er_opt(:, 2) == r);
    X_delim(N_parts) = m;
    Y_delim(N_parts) = C_opt(m);
    X_init_opt(N_parts) = X_init(r);
    Y_init_opt(N_parts) = Y_init(r);
end


diff = [];
m = 0;
for nd=X_delim
    if nd > 1 && ~isnan(C_opt(nd - 1))
        m = m+1;
        diff(m) = C_opt(nd) - C_opt(nd - 1);
    end
end

% Z = zeros(size(diff));
% figure;
% hold on;
% plot(diff);
% plot(Z);
% hold off;
% pause;

%% figures
figure;
imagesc((0:L-1)/L, (L/Nfft)*(1:Nfft), abs(TFR_it));
set(gca,'ydir','normal');
axis square
hold on;
plot((0:L-1)/L, C_opt*L/Nfft, 'r');
plot((X_delim-1)/L, Y_delim*L/Nfft, 'gx', 'MarkerSize', 10);
text((X_delim-1)/L, Y_delim*L/Nfft, string(E_parts), 'Color', 'w');
plot((X_init_opt-1)/L, Y_init_opt*L/Nfft, 'r+', 'MarkerSize', 10);
hold off;

%% display more info
% Cr(Cr == 0) = nan;
% figure;
% imagesc((0:L-1)/L, (L/Nfft)*(1:Nfft), abs(TFR));
% set(gca,'ydir','normal');
% axis square;
% hold on;
% plot((0:L-1)/L, Cr(:, :)*L/Nfft, 'r');
% hold off;

%% set to zero the ridge
% for n=1:L
%     B = real(q(c(n), n));
%     eta_k = ceil(Nfft/L*3/sqrt(2*pi)*sqrt(1/sigma^2 + sigma^2*B^2));
%     range = max(1, c(n)-eta_k):min(Nfft, c(n)+eta_k);
%     TFR_it(range, n) = 0;
% end

end