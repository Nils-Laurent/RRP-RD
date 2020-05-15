function [Cs, Qs, crossing, E_max] = novel_WPF_iterate(TFR, RRP_inter, RRP_E_inter, range_vec, Nr, degree_WPF)
Ni = length(RRP_E_inter);
[Nfft, L] = size(TFR);

Cs_it = zeros(Nr, L);
Qs_it = zeros(Nr, degree_WPF+1);
KY_lower_it = zeros(Nr, L);
KY_upper_it = zeros(Nr, L);

E_max = 0;
Cs = zeros(Nr, L);
Qs = zeros(Nr, degree_WPF+1);
crossing = 1;

ir_vec = zeros(1, Ni);
pr_vec = zeros(1, Ni);

SPEC = abs(TFR).^2;

[E_sorted, order_RS] = sort(RRP_E_inter, 'descend');
RRP_sorted = RRP_inter(order_RS, :);


R_inter = zeros(size(TFR));
for r=1:Ni
    for n=1:L
        kr = RRP_sorted(r, n);
        if kr > 0
            R_inter(kr, n) = E_sorted(r);
        end
    end
end

for ir=1:Ni
    % prepare initialization
    ir_vec(ir) = ir;
    [Cs_in_it, pr_new, N_init, valid_init] =...
        novel_WPF_init(TFR, RRP_sorted, E_sorted, ir_vec, pr_vec, Nr);
    pr_vec = pr_new;
    
    if N_init == 0
        if valid_init == 0
            ir_vec(ir) = 0;
        end
        continue;
    end
    
    % test polynomials
    E_new =  0;
    R_ip = R_inter;
    for p=1:Nr
        R_modes = zeros(size(TFR));
        for p2=1:Nr
            if p2 == p
                % keep RRP group concerning mode "p"
                continue;
            end
            
            % remove RRP group concerning "p2" neq "p" in input R_ip
            for n=1:L
                k = Cs_in_it(p2, n);
                if k > 0
                    R_modes(k, n) = R_ip(k, n);
                    R_ip(k, n) = 0;
                end
            end
        end
        
        [Cs_it(p, :), Qs_it(p, :), KY_lower_it(p, :), KY_upper_it(p, :), R_new] =...
            novel_WPF(R_ip, Cs_in_it(p, :), range_vec(p, :), degree_WPF);
        R_ip = R_new;
        R_ip = R_ip + R_modes;
    end
    
    %% check if modes are crossing in the TF plane
    crossing_it = 0;
    KY = zeros(size(KY_lower_it));
    for p=1:Nr
        KY(p, :) = polyval(Qs_it(p, :), (0:(L-1))/L);
        KY(p, :) = round(KY(p, :)*Nfft/L) + 1;
        KY(p, :) = min(Nfft, max(1, KY(p, :)));
    end
    for p=1:Nr
        cmp = KY(p, :);
        cmp(cmp == 1) = 0;
        cmp(cmp == Nfft) = Nfft - 1;
        for p2=(p+1):Nr
            if sum(cmp >= KY(p2, :)) > 0
                crossing_it = crossing_it + 1;
            end
        end
    end
    
    %% compute energy with proper intervals
    KYd = KY_lower_it;
    KYu = KY_upper_it;
    for n=1:L
        for p=1:Nr
            % if the size is of I_LC is too small, continue
            if KYd(p, n) == Nfft || KYu(p, n) == 1
                continue;
            end
            
            ensure_v = 1;
            % skip the case HIGH FREQ region of below LOW FREQ
            if p == Nr || KYu(p + 1, n) <= KYd(p, n)
                ensure_v = 0;
            end
            
            % ensure disjoint intervals
            if ensure_v == 1 && KYu(p, n) >= KYd(p+1, n)
                mid = KYd(p+1, n) + round((KYu(p, n) - KYd(p+1, n))/2);
                KYu(p, n) = mid;
                KYd(p+1, n) = mid+1;
            end
            
            E_new = E_new + sum(SPEC(KYd(p, n):KYu(p, n), n));
        end
    end
    
    %% register max info
    if E_new > E_max && crossing_it <= crossing
        E_max = E_new;
        Qs = Qs_it;
        Cs = Cs_it;
        
        crossing = crossing_it;
    end
end

fprintf('crossing = %u\n', crossing);

end

