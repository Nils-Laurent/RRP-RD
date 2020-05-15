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

    pr_new = pr_in;
    for n=1:L
        p_id = 0;
        k_pr = zeros(1, Nr);
        r_pr = zeros(1, Nr);
        for r=ir_vec(ir_vec > 0)
            k = RRP_sorted(r, n);
            if k > 0
                p_id = p_id + 1;
                if p_id > Nr
                    break;
                end
                k_pr(p_id) = k;
                r_pr(p_id) = r;
            end
        end

        if p_id == Nr
            [~, p_id_sort] = sort(k_pr);
            for p=1:Nr
                pp_id = p_id_sort(p);
                r = r_pr(pp_id);
                if pr_new(r) == 0 || pr_new(r) == p
                    pr_new(r) = p;
                else
                    mode_conflict = 1;
                end
            end
        end

        if mode_conflict == 1 || p_id > Nr
            break;
        end

        if p_id == Nr
            N_init = N_init + 1;
        end
    end
    
    if mode_conflict == 1
        valid_init = 0;
%         fprintf("init: multiple modes on RRP\n");
    end
    
    if p_id > Nr
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
%     
%     % interpolation energy matrix
%     for r=1:ir
%         for n=1:L
%             kr = RRP_sorted(r, n);
%             if kr > 0
%                 R_ir(kr, n) = E_sorted(r);
%             end
%         end
%     end
end

