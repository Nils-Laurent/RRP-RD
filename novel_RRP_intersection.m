function [RRP_merged, RRP_E_merged, RRP_M_ES] = novel_RRP_intersection(TFR, RRP_stable, RRP_E_stable, range_I, LNh)

% fprintf("merge related structures\n");

[Nfft, L] = size(TFR);
Ns = length(RRP_E_stable);
r_rel_vec = zeros(1, Ns);
Nm = 0;

function [TFR_r, min_n, max_n, min_k, max_k] = RRP_supp(Nfft, L, k0_vec, LNh, range_I)
    nf = 1;
    
    TFR_r = [];
    min_k = Nfft;
    max_k = 1;
    
    while k0_vec(nf) == 0
        if nf == L
            return;
        end
        nf = nf + 1;
    end
    min_n = max(1, nf - LNh);
    
    if nf - 1 > 0
        nA = max(1, nf - LNh);
        nB = max(1, nf - 1);
        kA = max(1, k0_vec(nf) - range_I);
        kB = min(Nfft, k0_vec(nf) + range_I);
        A = zeros(Nfft,nB-nA+1); 
        A(kA:kB,:) = 1;
        TFR_r = A; 
    end
    
    while k0_vec(nf) > 0
        kA = max(1, k0_vec(nf) - range_I);
        kB = min(Nfft, k0_vec(nf) + range_I);
        
        min_k = min(min_k, kA);
        max_k = max(max_k, kB);
        
        A = zeros(Nfft,1); 
        A(kA:kB) = 1;
        
        if isempty(TFR_r)
            TFR_r = A;
        else
            TFR_r = [TFR_r A];
        end
        
        if nf == L
            max_n = L;
            TFR_r = TFR_r(min_k:max_k,:);
            return;
        end
        nf = nf + 1;
    end
    max_n = min(L, nf - 1 + LNh);
    
    if nf <= L
        nA = min(L, nf);
        nB = min(L, nf + LNh - 1);
        kA = max(1, k0_vec(nf-1) - range_I);
        kB = min(Nfft, k0_vec(nf-1) + range_I);
        A = zeros(Nfft,nB-nA+1); 
        A(kA:kB,:) = 1;
        TFR_r = [TFR_r A];
    end
   TFR_r = TFR_r(min_k:max_k,:);
end

[~, rE_vec] = sort(RRP_E_stable, 'descend');

for r0=rE_vec
    
    if r_rel_vec(r0) > 0
        continue;
    end
    
    r_rel_vec(r0) = r0;
    
    k0_vec = RRP_stable(r0, :);
    supp_r0 = zeros(Nfft, L);
    [A_r0, n1, n2, k1, k2] = RRP_supp(Nfft, L, k0_vec, LNh, range_I);
    supp_r0(k1:k2, n1:n2) = A_r0;
    
    r_vec = zeros(1, Ns);
    
    %% loop on other RRP
    isec = 1;
    while isec == 1
        isec = 0;
        for r=1:Ns
            if r_rel_vec(r) > 0 || r_vec(r) > 0
                continue;
            end
            [supp_r, n1, n2, k1, k2] = RRP_supp(Nfft, L, RRP_stable(r, :), LNh, range_I);
            
            supp_I = supp_r0(k1:k2, n1:n2).*supp_r;
            if sum(supp_I(:)) > 0
                r_vec(r) = r;
                supp_r0(k1:k2, n1:n2) = max(supp_r0(k1:k2, n1:n2), supp_r);
                isec = 1;
                break;
            end
            
        end
    end
    
%     figure;
%     imagesc(1:L, 1:Nfft, supp_r0);
%     set(gca,'ydir','normal');
%     colormap(flipud(gray));
%     axis square
%     title("supp_r0");
%     pause;
    
    %% merge process
%     fprintf("intersection : r0 = %u\n", r0);
    r_vec = r_vec(r_vec > 0);
    for r=r_vec
        r_rel_vec(r) = r0;
    end
    
    Nm = Nm + 1;
end

RRP_merged = zeros(Nm, L);
RRP_E_merged = zeros(1, Nm);
RRP_M_ES = zeros(1, Nm);

rm = 1;
for r0=1:Ns
    if r0 ~= r_rel_vec(r0)
        continue;
    end
    
    RRP_M_ES(rm) = RRP_E_stable(r0);
    
    e_vec = zeros(1, Ns);
    for r=1:Ns
        if r_rel_vec(r) == r0
            e_vec(r) = RRP_E_stable(r);
        end
    end
    RRP_E_merged(rm) = sum(e_vec);
    
    [~, r_sorted] = sort(e_vec, 'descend');
    for n=1:L
        for rs=r_sorted
            if e_vec(rs) == 0
                break;
            end
            if RRP_stable(rs, n) > 0
                RRP_merged(rm, n) = RRP_stable(rs, n);
                break;
            end
        end
    end
    
    rm = rm + 1;
end

%% plot merge TFR

% TFR_merged = zeros(size(TFR));
% for r=1:Nm
%     for n=1:L
%         kr = RRP_merged(r, n);
%         if kr > 0
%             TFR_merged(kr, n) = RRP_E_merged(r);
%         end
%     end
% end
% 
% figure;
% imagesc(1:L, 1:Nfft, TFR_merged);
% set(gca,'ydir','normal');
% axis square
% colormap(flipud(gray));
% title("TFR intersection");

end

