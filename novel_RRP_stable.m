function [RRP_stable, RRP_E_stable] = novel_RRP_stable(TFR, S_LM, QM, sigma_s, Nr, scale)
fprintf('compute ridge portions\n');
scale = 2*ceil(scale/2);

[Nfft, L] = size(S_LM);

Sr = zeros([Nr, L, L]);

fprintf('scale = %u :', scale);
for m=1:L
    if mod(m, 512) == 0
        fprintf(' %d', m);
    end
    
    [~, k_vec] = sort(S_LM(:, m), 'descend');
    k_vec_Nr = k_vec(1:Nr);
    k_vec_Nr = sort(k_vec_Nr);
    for p=1:Nr
        kp = k_vec_Nr(p);
        if S_LM(kp, m) == 0
            break;
        end
        
        [Sr(p, m, :), ~] = novel_partial_RD(S_LM, m, kp, QM, Nr, sigma_s);
    end
end
fprintf('\n');

% DISP = zeros(L, L);
% DISP(:, :) = Sr(1, :, :);
% figure;
% imagesc(1:L, 1:L, DISP);
% set(gca,'ydir','normal');
% axis square
% colormap(flipud(gray));
% title("Sr m1");

%% compute stability
% fprintf("compute stability\n");

S_stable = Sr;
Stable_m0_vec = zeros(Nr, L);
Ns = 0;
for p=1:Nr
    m0 = 1;
    V_ref = S_stable(p, m0, :);
    C = 1;
    for m=2:L
        if sum(V_ref == S_stable(p, m, :)) == L
            C = C + 1;
        else
            if C <= scale
                for q=m0:(m-1)
                    S_stable(p, q, :) = zeros(1, L);
                end
            else
                Stable_m0_vec(p, m0) = 1;
            end
            m0 = m;
            V_ref = S_stable(p, m0, :);
            C = 1;
            Ns = Ns + 1;
        end
    end
end

%% create stable RRP structure
% fprintf("create light stable structure\n");
RRP_stable = zeros(Ns, L);
RRP_m0 = zeros(1, Ns);
r = 1;
for p=1:Nr
    for m=1:L
        if Stable_m0_vec(p, m) == 0
            continue;
        end
        
        RRP_stable(r, :) = S_stable(p, m, :);
        %% do not allow zero RRP
        if sum(RRP_stable(r, :) == zeros(1, L)) == L
            continue;
        end
        
        RRP_m0(r) = m;
        
        %% do not allow duplicates
        for r_test=1:(r-1)
            if sum(RRP_stable(r_test, :) == RRP_stable(r, :)) == L
                RRP_m0(r) = 0;
                break;
            end
        end
        
        if RRP_m0(r) > 0
            Ns = r;
            r = r + 1;
        end
    end
end
RRP_stable = RRP_stable(1:Ns, :);
% RRP_m0 = RRP_m0(1:Ns);

RRP_E_stable = zeros(1, Ns);
for r=1:Ns
    kr_vec = RRP_stable(r, :);
    for n=1:L
        if kr_vec(n) == 0
            continue;
        end
        RRP_E_stable(r) = RRP_E_stable(r) + abs(TFR(kr_vec(n), n))^2;
    end
end

%% show energy TFR
% R_Energy = zeros(Nfft, L);
% for r=1:Ns
%     kr_vec = RRP_stable(r, :);
%     for n=1:L
%         if kr_vec(n) == 0
%             continue;
%         end
%         R_Energy(kr_vec(n), n) = RRP_E_stable(r);
%     end
% end
% 
% figure;
% imagesc(1:L, 1:Nfft, R_Energy);
% set(gca,'ydir','normal');
% axis square
% colormap(flipud(gray));
% title("R energy");

end

