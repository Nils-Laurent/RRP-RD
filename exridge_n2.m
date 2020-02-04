function [C_opt] = exridge_new(TFR, Lg, sigma, q, omega, omega2, C)

[Nfft, L] = size(TFR);

gamma = median(abs(real(TFR(:))))/0.6745;

% Values from the Chi-squared table : P[X > Cp] = p
%ratio = 1/100;
%C1 = 9.2103;
ratio = 1/10;
C10 = 4.6052;

TH = gamma*sqrt(C10);

SZs = [floor(Lg/2), Lg, Lg+floor(Lg/2), 2*Lg+1, 3*Lg, 4*Lg];
N_SZ = length(SZs);
Wmat = zeros(N_SZ, L);
Cmat = zeros(N_SZ, L);

N_shift = 8;
shift_sz = max(2, floor((2*Lg+1)/N_shift));
Shifts = (0:shift_sz:8*shift_sz) - 4*shift_sz;

for step_sz=SZs
    Steps = 1:step_sz:L;

    for step=Steps
        for shift=Shifts
            n0 = max(1, step+shift);
            n0 = min(n0, L);
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
    end
    
    Wmat(step_sz, :);
    Cmat(step_sz, :);
end

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

%% Look for energy where response is missing
C_lin_FM = C_opt;
C_lin = C_opt;
C_lin_init = C_opt;

C1it_lin = C_opt;
C2it_lin = C_opt;
R = max(2, floor(L/(4*Nfft)));

Y = abs(TFR);
for m=1:M
    n1 = X1_opt(m) +1;
    n2 = X2_opt(m) -1;
    
    %% high response extremities fitting
    nf1 = n1;
    nf2 = n2;
    line_E = 0;
    line_mi = round(linspace(C_opt(n1-1), C_opt(n2+1), n2-n1+1));
    line_m = line_mi;
    for n1t=max(2, n1-Lg):n1
        if isnan(C_opt(n1t-1))
            continue;
        end
        for n2t=n2:min(L-1, n2+Lg)
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

end