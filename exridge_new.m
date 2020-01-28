function [Cs_lin] = exridge_new(TFR, Lg, sigma, q, C)

[Nfft, L] = size(TFR);
%c = zeros(L, 1);

gamma = median(abs(real(TFR(:))))/0.6745;

% Values from the Chi-squared table : P[X > Cp] = p
%C1 = 9.2103;
ratio = 1/10;
C10 = 4.6052;

TH = gamma*sqrt(C10);

%N_lin = min(60,floor(L/8));
%init_set = floor(linspace(L/(N_lin+1), L-L/(N_lin+1), N_lin));
%init_set = floor(linspace(L/(N_lin+1), L-L/(N_lin+1), N_lin));
N_lin = floor(L/(2*Lg));
step = floor(L/N_lin);
init_set = 1:step:(L-step);
N_lin = length(init_set);
%R = floor(min(step, Lg/16));
R = max(2, floor(L/(4*Nfft)));
N_shift = floor(step/R);
shifts = 0:R:(N_shift*R-1);
N_shift = length(shifts)


Cs_shift = zeros(N_shift, L);
Es_shift = zeros(N_shift);
%Er = zeros(N_lin, N_shift);

for n_shift = 1:N_shift
    Cs_lin = zeros(N_lin, L);
    Es_lin = zeros(N_lin, 1);
    E_matrix = zeros(N_lin, L);
    X_init = zeros(N_lin, 1);
    Y_init = zeros(N_lin, 1);
    Ar = zeros(N_lin, 1);
    Br = zeros(N_lin, 1);
    for n_lin = 1:N_lin
        n0 = init_set(n_lin) + shifts(n_shift);
        [e0, y0] = max(abs(TFR(:, n0)).*(abs(TFR(:, n0))>3*gamma));
        if e0 == 0
            continue;
        end
            
        X_init(n_lin) = n0;
        %Y_init(n_lin, n_shift) = y0;
        E_matrix(n_lin, n0) = e0^2;
        e = e0^2;
        Cs_lin(n_lin, n0) = y0;
        Ar(n_lin) = n0;
        Br(n_lin) = n0;

        RQ = zeros(L, 1);
        RQ(n0) = round(Nfft/(L^2)*real(q(Cs_lin(n_lin, n0), n0)));

        %% forward iterations
        for n=(n0+1):L
            next_index = min(Nfft, Cs_lin(n_lin, n-1) +RQ(n-1));
            next_index = max(1, next_index);
            Im = max(1, next_index -C):min(Nfft, next_index +C);

            Im_ratio = sum(abs(TFR(Im, n)) > TH)/length(Im);
            if Im_ratio <= ratio
                break;
            end

            [v_arg, arg] = max(abs(TFR(Im, n)));
            E_matrix(n_lin, n) = v_arg^2;
            e = e + v_arg^2;
            Br(n_lin) = n;
            Cs_lin(n_lin, n) = arg + Im(1)-1;
            RQ(n) = round(Nfft/(L^2)*real(q(Cs_lin(n_lin, n), n)));
        end

        %% backward iterations
        for n=(n0-1):-1:1
            next_index = min(Nfft, Cs_lin(n_lin, n+1) -RQ(n+1));
            next_index = max(1, next_index);
            Im = max(1, next_index -C):min(Nfft, next_index +C);

            Im_ratio = sum(abs(TFR(Im, n)) > TH)/length(Im);
            if Im_ratio <= ratio
                break;
            end

            [v_arg, arg] = max(abs(TFR(Im, n)));
            E_matrix(n_lin, n) = v_arg^2;
            e = e + v_arg^2;
            Ar(n_lin) = n;
            Cs_lin(n_lin, n) = arg + Im(1)-1;
            RQ(n) = round(Nfft/(L^2)*real(q(Cs_lin(n_lin, n), n)));
        end

        %Er(n_lin, n_shift) = e;
        Es_lin(n_lin) = e;
    end
    
%     [Es_lin, Er_order] = sort(Es_lin, 'descend');
%     Cs_lin = Cs_lin(Er_order, :);
%     X_init = X_init(Er_order);
%     Y_init = Y_init(Er_order);
%     Ar = Ar(Er_order);
%     Br = Br(Er_order);
%     
%     %% compute C_opt(n_off)
%     opt_map = zeros(L, 1);
%     C_opt = zeros(L, 1);
%     Nr_opt = 0;
%     for n_lin=1:N_lin
%         nd = X_init(n_lin);
%         if opt_map(nd) > 0
%             continue;
%         end
% 
%         Nr_opt = Nr_opt + 1;
% 
%         for m = Ar(n_lin):Br(n_lin)
%             if opt_map(m) == 0
%                 opt_map(m) = n_lin;
%                 C_opt(m) = Cs_lin(opt_map(m), m);
%             end
%             if C_opt(m) > 0
%                 Es_shift(n_shift) = Es_shift(n_shift) + abs(TFR(C_opt(m), m))^2;
%             end
%         end
%     end

    C_opt_shift = zeros(L, 1);
    
    for n=1:L
        [v, arg_lin] = max(E_matrix(:, n));
        if v > 0
            C_opt_shift(n) = Cs_lin(arg_lin, n);
        else
            C_opt_shift(n) = 0;
        end
    end
    
    Cs_shift(n_shift, :) = C_opt_shift;
    
%     figure;
%     imagesc((0:L-1)/L, (L/Nfft)*(1:Nfft), abs(TFR_it));
%     set(gca,'ydir','normal');
%     axis square
%     hold on;
%     plot((0:L-1)/L, (L/Nfft)*(C_opt -  1), 'r');
%     pause
end

%% compare offset dependent ridges
%Cs_test = zeros(N_shift, L);
% test_std = zeros(L, 1);
% C_opt_shift = zeros(L, 1);
% for n=1:L
%     A = 0;
%     for n_shift=1:N_shift
%         A = A + Cs_shift(n_shift, n);
%         if isnan(A)
%             break;
%         end
%     end
% %     X = nonzeros(Cs_shift(:,n));
% %     if isempty(X)
% %         continue;
% %     end
% %     X2 = unique(X);
% %     ns = histc(X, X2);
% %     [~, arg] = max(ns);
% %     y = X2(arg);
%     cur_std = std(Cs_shift(:, n));
%     test_std(n) = cur_std;
%     TH_std = sqrt((Nfft^2-1)/12);
%     if (cur_std - TH_std) > 0 || length(X) < 10
%         y = nan;
%     else
%         y = median(X);
%     end
%     if length(X) >= 10,
%      y = median(X);
%     else
%      y = nan;
%     end 
    %if n == 1775
    % Cs_shift(:,n)
    % pause
    %end 
    %[~, arg] = max(abs(TFR(Cs_shift(:, n), n)));
%     C_opt_shift(n) = y;
% end

C_opt = zeros(L, 1);
max_count = zeros(L, 1);
len_X = zeros(L, 1);
for n=1:L
    X = nonzeros(Cs_shift(:,n));
    len_X(n) = length(X);
    if isempty(X)
        continue;
    end
    X2 = unique(X);
    [v, arg] = max(histc(X, X2));
    % freq cos : v > 15
    %if v > floor(2*N_shift/3)
    if v > 20
        max_count(n) = v;
        y = X2(arg);
        C_opt(n) = y;
    else
        C_opt(n) = nan;
    end
end

C_opt(C_opt == 0) = nan;

% figure;
% imagesc(0:L-1, (1:Nfft), abs(TFR));
% set(gca,'ydir','normal');
% axis square
% hold on;
% plot(0:L-1, (Cs_shift(:, :)));
% hold off;
% pause

figure;
imagesc(0:L-1, (1:Nfft), abs(TFR));
set(gca,'ydir','normal');
axis square
hold on;
plot(0:L-1, C_opt, 'r');
hold off;

%% look for ridge parts
S = 0;
for n=1:L
    if ~isnan(C_opt(n))
        S = n;
        break;
    end
end

X1_opt = zeros(L, 1);
X2_opt = zeros(L, 1);
m = 1;
M = 0;
prev = C_opt(S);
for n = (S+1):L
    cur = C_opt(n);
    if ~isnan(prev) && isnan(cur)
        X1_opt(m) = (n-1);
    elseif isnan(prev) && ~isnan(cur)
        X2_opt(m) = n;
        M = m;
        m = m + 1;
    end
    prev = cur;
end
X1_opt = nonzeros(X1_opt);
X2_opt = nonzeros(X2_opt);

%% ridge completion
for m=1:M
    n1 = X1_opt(m);
    n2 = X2_opt(m);
    k1 = C_opt(n1);
    k2 = C_opt(n2);
    ka = min(k1, k2);
    kb = max(k1, k2);
    if kb - ka < 2
        for n=(n1+1):(n2-1)
            C_opt(n) = kb;
        end
        continue;
    end
    
    Y = abs(TFR(ka:kb, (n1+1):(n2-1)));
    Y_max_loc = zeros(1,n2-n1-1);
    
%     for n=2:(n2-n1)
%         for k=2:(kb-ka)
%             if Y(k-1, n) < Y(k, n) && Y(k+1, n) < Y(k, n)
%                 Y_max_loc(k, n) = Y(k, n);
%             end
%         end
%     end
%     Y = nonzeros(Y_max_loc(:));
    for n=1:n2-n1-1
     Y_max_loc(n) = max(Y(:,n));
    end
    
    %Y_max_loc = Y(Y
    Y2 = unique(Y(:));
    [~, arg] = max(histc(Y(:), Y2));
    ridge_v = Y2(arg);
    for n=(n1+1):(n2-1)
        %[~, arg] = min(abs(abs(TFR(ka:kb, n)) - ridge_v));
        [~, arg] = max(abs(TFR(ka:kb, n)));
        C_opt(n) = ka + arg -1;
    end
end

figure;
imagesc(0:L-1, (1:Nfft), abs(TFR));
set(gca,'ydir','normal');
axis square
hold on;
plot(0:L-1, C_opt, 'r');
hold off;

%% figures
% figure;
% imagesc((0:L-1)/L, (L/Nfft)*(1:Nfft), abs(TFR_it));
% set(gca,'ydir','normal');
% axis square
% hold on;
% plot((X_init - 1)/L, (L/Nfft)*(Y_init -  1));
% plot((0:L-1)/L, C_opt*L/Nfft, 'r');
% plot((X_delim-1)/L, Y_delim*L/Nfft, 'gx', 'MarkerSize', 10);
% text((X_delim-1)/L, Y_delim*L/Nfft, string(E_parts), 'Color', 'w');
% plot((X_init_opt-1)/L, Y_init_opt*L/Nfft, 'r+', 'MarkerSize', 10);
% hold off;

end