function [Cr_in, pr_new, N_init, valid_init] =...
    novel_WPF_init(TFR, RRP_sorted, E_sorted, ir_vec, pr_in, Nr)
    [~, L] = size(TFR);
%     R_ir = zeros(size(TFR));
    Cr_in = zeros(Nr, L);
    
    Nm = length(ir_vec);
    ir = max(ir_vec);
    
    valid_init = 1;
    N_init = 0;
    mode_conflict = 0;

    % assign frequency to each mode when Nr responses exist
    % define (frequency, RRP) couple (k, r) at each n in [1, ..., L]
    pr_new = pr_in;
    for n=1:L
        p_id = 0; % identifier for one of the modes p in [1, ..., Nr]
        k_pr = zeros(1, Nr);
        r_pr = zeros(1, Nr);
        for r=ir_vec(ir_vec > 0)
            k = RRP_sorted(r, n);
            if k > 0
                p_id = p_id + 1;
                if p_id > Nr
                    % too many responses at time index n
                    break;
                end
                
                % define couple (k, r)
                k_pr(p_id) = k;
                r_pr(p_id) = r;
            end
        end
        
        % number of responses at time index n
        N_resp = p_id;

        if N_resp == Nr
            % there exist one RRP for each mode
            % assign RRP to modes from LF to HF
            [~, p_id_sort] = sort(k_pr);
            for p=1:Nr
                % access mode information from LF to HF
                p_id = p_id_sort(p);
                r = r_pr(p_id);
                
                if pr_new(r) == 0 || pr_new(r) == p
                    % assign RRP with TF (n, k) to mode p
                    pr_new(r) = p;
                else
                    % RRP r was already assigned => conflict
                    mode_conflict = 1;
                end
            end
        end

        if mode_conflict == 1 || N_resp > Nr
            break;
        end

        if N_resp == Nr
            N_init = N_init + 1;
        end
    end
    
    if mode_conflict == 1
        valid_init = 0;
%         fprintf("init: multiple modes on RRP\n");
    end
    
    if N_resp > Nr
        valid_init = 0;
%         fprintf("init: too many responses\n");
    end
    
    if valid_init == 0
        N_init = 0;
        pr_new = pr_in;
        return;
    end
    
    % assign Cr_in
    for r=1:Nm
        p = pr_new(r);
        if p > 0
            Cr_in(p, :) = max(Cr_in(p, :), RRP_sorted(r, :));
        end
    end
    
    if N_init == 0
        return;
    end
    
    % interpolation energy matrix
%     R_ir = zeros(size(TFR));
%     for r=1:ir
%         for n=1:L
%             kr = RRP_sorted(r, n);
%             if kr > 0
%                 R_ir(kr, n) = E_sorted(r);
%             end
%         end
%     end
%     
%     [Nfft, ~] = size(TFR);
%     figure;
%     imagesc(1:L, 1:Nfft, R_ir);
%     set(gca,'ydir','normal');
%     axis square
%     colormap(flipud(gray));
%     title("R ir");
%     pause;
end

