function [RRP_merged, RRP_E_merged, RRP_M_ES] = novel_RRP_intersection(TFR, RRP_stable, RRP_E_stable, range_I, LNh)

% fprintf("merge related structures\n");

[Nfft, L] = size(TFR);
Ns = length(RRP_E_stable);
r_rel_vec = zeros(1, Ns);
Nm = 0;

function [TFR_r] = RRP_supp(Nfft, L, k0_vec, LNh, range_I)
    TFR_r = zeros(Nfft, L);
    nf = 1;
    
    while k0_vec(nf) == 0
        if nf == L
            return;
        end
        nf = nf + 1;
    end
    
    if nf - 1 > 0
        nA = max(1, nf - LNh);
        nB = max(1, nf - 1);
        kA = max(1, k0_vec(nf) - range_I);
        kB = min(Nfft, k0_vec(nf) + range_I);
        TFR_r(kA:kB, nA:nB) = 1;
    end
    
    while k0_vec(nf) > 0
        kA = max(1, k0_vec(nf) - range_I);
        kB = min(Nfft, k0_vec(nf) + range_I);
        TFR_r(kA:kB, nf) = 1;
        
        if nf == L
            return;
        end
        nf = nf + 1;
    end
    
    if nf <= L
        nA = min(L, nf);
        nB = min(L, nf + LNh);
        kA = max(1, k0_vec(nf-1) - range_I);
        kB = min(Nfft, k0_vec(nf-1) + range_I);
        TFR_r(kA:kB, nA:nB) = 1;
    end
end

[~, rE_vec] = sort(RRP_E_stable, 'descend');

for r0=rE_vec
    
    if r_rel_vec(r0) > 0
        continue;
    end
    
    r_rel_vec(r0) = r0;
    
    k0_vec = RRP_stable(r0, :);
    supp_r0 = RRP_supp(Nfft, L, k0_vec, LNh, range_I);
    
    r_vec = zeros(1, Ns);
    
    %% loop on other RRP
    isec = 1;
    while isec == 1
        isec = 0;
        for r=1:Ns
            if r_rel_vec(r) > 0 || r_vec(r) > 0
                continue;
            end
            supp_r = RRP_supp(Nfft, L, RRP_stable(r, :), LNh, range_I);
            
            supp_I = supp_r0.*supp_r;
            if sum(supp_I(:)) > 0
                r_vec(r) = r;
                supp_r0 = max(supp_r0, supp_r);
                isec = 1;
                break;
            end
        end
    end
    
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

