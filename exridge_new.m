function [C_opt] = exridge_new(TFR, Lg, sigma, q, omega, omega2, C)

[Nfft, L] = size(TFR);

gamma = median(abs(real(TFR(:))))/0.6745;

% Values from the Chi-squared table : P[X > Cp] = p
%ratio = 1/100;
%C1 = 9.2103;
ratio = 1/10;
C10 = 4.6052;

TH = gamma*sqrt(C10);

%N_lin = min(60,floor(L/8));
%init_set = floor(linspace(L/(N_lin+1), L-L/(N_lin+1), N_lin));
N_lin = floor(L/(2*Lg));
step = floor(L/N_lin);
init_set = 1:step:(L-step);
N_lin = length(init_set);
R = max(2, floor(L/(4*Nfft)));
N_shift = floor(step/R);
shifts = 0:R:(N_shift*R-1);
N_shift = length(shifts);


Cs_shift = zeros(N_shift, L);

for n_shift = 1:N_shift
    Cs_lin = zeros(N_lin, L);
    Es_lin = zeros(N_lin, 1);
    E_matrix = zeros(N_lin, L);
    X_init = zeros(N_lin, 1);
    Ar = zeros(N_lin, 1);
    Br = zeros(N_lin, 1);
    for n_lin = 1:N_lin
        n0 = init_set(n_lin) + shifts(n_shift);
        [e0, y0] = max(abs(TFR(:, n0)).*(abs(TFR(:, n0))>3*gamma));
        if e0 == 0
            continue;
        end
            
        X_init(n_lin) = n0;
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

        Es_lin(n_lin) = e;
    end

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

%% detect high response from Cs_shift
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
    % if v > 15
    % if v > 20
    if v > floor(2*N_shift/3)
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
% plot(0:L-1, C_opt, 'r');
% hold off;

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

%% compute pchip from C_opt
% Xpchip = zeros(L, 1);
% Ypchip = zeros(L, 1);
% for n=1:L
%     if ~isnan(C_opt(n))
%         Xpchip(n) = n;
%         Ypchip(n) = C_opt(n);
%     end
% end
% Xpchip = nonzeros(Xpchip);
% Ypchip = nonzeros(Ypchip);
% C_pchip = pchip(Xpchip, Ypchip, 1:L);

%% C_opt_pchip
% C_opt_pchip = C_opt;
% Y = abs(TFR);
% for m=1:M
%     n1 = X1_opt(m) +1;
%     n2 = X2_opt(m) -1;
%     for n=n1:n2
%         k_pc = C_pchip(n);
%         k_max = inf;
%         for k=2:(Nfft-1)
%             if Y(k-1, n) < Y(k, n) && Y(k+1, n) < Y(k, n)
%                 if abs(k - k_pc) < abs(k_max - k_pc)
%                     k_max = k;
%                 end
%             end
%         end
%         C_opt_pchip(n) = k_max;
%     end
% end
% 
% figure;
% imagesc((0:L-1)/L, (1:Nfft), abs(TFR));
% set(gca,'ydir','normal');
% axis square
% hold on;
% plot((0:L-1)/L, C_pchip, 'k');
% plot((0:L-1)/L, C_opt_pchip, 'r');
% hold off;

%% C_opt_lin
% C_opt_lin_w2 = C_opt;
% C_opt_lin_w = C_opt;
%C_opt_lin = C_opt;
C1it_lin = C_opt;
C2it_lin = C_opt;
C_lin_FM = C_opt;
C_lin_w = C_opt;
C_lin_w2 = C_opt;
C_lin = C_opt;
C_lin_init = C_opt;
Y = abs(TFR);
R = max(2, floor(L/(4*Nfft)));
for m=1:M
    n1 = X1_opt(m) +1;
    n2 = X2_opt(m) -1;
    
    %% high response extremities fitting
    nf1 = n1;
    nf2 = n2;
    line_E = 0;
    line_mi = round(linspace(C_opt(n1-1), C_opt(n2+1), n2-n1+1));
    line_m = line_mi;
    %if ~isnan(C_opt(n1-Lg-1)) && ~isnan(C_opt(n2+Lg+1))
    for n1t=n1-Lg:n1
        if isnan(C_opt(n1t-1))
            continue;
        end
        for n2t=n2:n2+Lg
            if isnan(C_opt(n2t+1))
                continue;
            end
            line_t = round(linspace(C_opt(n1t-1), C_opt(n2t+1), n2t-n1t+1));
            line_tE = 0;
            for nt=n1:n2
                line_tE = line_tE + Y(line_t(nt -n1t +1), nt)^2;
            end
            if line_tE > line_E
                nf1 = n1t;
                nf2 = n2t;
                line_E = line_tE;
                line_m = line_t;
            end
        end
    end

%     figure;
%     imagesc(1:L, 1:Nfft, abs(TFR));
%     set(gca,'ydir','normal');
%     axis square
%     hold on;
%     plot([ni1-1, ni2+1], [C_opt(ni1-1), C_opt(ni2+1)], 'r');
%     plot([n1-1, n2+1], [C_opt(n1-1), C_opt(n2+1)], 'g--');
%     hold off;
%     legend;
%     pause

%     fprintf("m = %f -------------------\n", m);
%     fprintf("n1 = %f -------------------\n", n1);
%     fprintf("diff_n1 = %f\n", n1-nf1);
%     fprintf("diff_n2 = %f\n", nf2-n2);

    line_m = line_m(n1-nf1+1:n2-nf1+1);
    C_lin(n1:n2) = line_m;
    C_lin_init(n1:n2) = line_mi;
    
    %% direct method
    
    C_lin_FM(n1:n2) = line_m;
    C_lin_w(n1:n2) = line_m;
    C_lin_w2(n1:n2) = line_m;
    for n=n1:n2
        k_line = line_m(n-n1+1);

        kw2 = round(omega2(k_line, n)*(Nfft/L))+1;
        kw = round(omega(k_line, n)*(Nfft/L))+1;

        [~, II] = min([abs(kw2 - k_line) abs(kw - k_line)]);
        if II == 1
            Cn = kw2;
        else
            Cn = kw;
        end
        C_lin_FM(n) = Cn;
        C_lin_w(n) = kw;
        C_lin_w2(n) = kw2;
    end

%     %% iteration init.
%     C1it_lin(n1:n2) = line_m;
%     C2it_lin(n1:n2) = line_m;
%     
%     %% iteration : fixed starting point
%     for p = n2:-R:n1
%         line_mp = round(linspace(C1it_lin(n1-1), C1it_lin(p+1), p-n1+1));
%         for n=n1:p
%             k_line = line_mp(n-n1+1);
% 
%             kw2 = round(omega2(k_line, n)*(Nfft/L))+1;
%             kw = round(omega(k_line, n)*(Nfft/L))+1;
% 
%             [~, II] = min([abs(kw2 - k_line) abs(kw - k_line)]);
%             if II == 1
%                 Cn = kw2;
%             else
%                 Cn = kw;
%             end
%             
%             if isnan(C1it_lin(n)) || Y(Cn) > Y(C1it_lin(n))
%                 C1it_lin(n) = Cn;
%             end
%         end
%     end
% 
%     %% iteration : fixed ending point
%     for p = n1:R:n2
%         line_mp = round(linspace(C2it_lin(p-1), C2it_lin(n2+1), n2-p+1));
%         for n=p:n2
%             k_line = line_mp(n-p+1);
% 
%             kw2 = round(omega2(k_line, n)*(Nfft/L))+1;
%             kw = round(omega(k_line, n)*(Nfft/L))+1;
% 
%             [~, II] = min([abs(kw2 - k_line) abs(kw - k_line)]);
%             if II == 1
%                 Cn = kw2;
%             else
%                 Cn = kw;
%             end
%             
%             if isnan(C2it_lin(n)) || Y(Cn) > Y(C2it_lin(n))
%                 C2it_lin(n) = Cn;
%             end
%         end
%     end
end

figure;
imagesc(1:L, 1:Nfft, abs(TFR));
set(gca,'ydir','normal');
axis square
title("initialization");
hold on;
plot(1:L, C_lin, 'r');
hold off;

figure;
imagesc(1:L, 1:Nfft, abs(TFR));
set(gca,'ydir','normal');
axis square
title("positionnement direct (min)");
hold on;
plot(1:L, C_lin, 'k');
plot(1:L, C_lin_FM, 'r--');
hold off;

% figure;
% imagesc(1:L, 1:Nfft, abs(TFR));
% set(gca,'ydir','normal');
% axis square
% title("positionnement iteratif");
% hold on;
% plot(1:L, C_lin, 'k');
% plot(1:L, C1it_lin, 'r--');
% plot(1:L, C2it_lin, 'g--');
% hold off;

C_opt = C_lin_FM;

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