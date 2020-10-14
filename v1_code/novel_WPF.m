function [C_new, p, KY_lower, KY_upper, R_out] = novel_WPF(R_energy, Cr_init, range_vec, degree_WPF)

[Nfft, L] = size(R_energy);

%% polyfit RD
C_prev = zeros(1, L);
C_new = Cr_init;

ITER = 0;
while sum(C_prev == C_new) < L
    ITER = ITER + 1;
    m = 0;
    iPX = [];
    iPY = [];
    iPW = [];
    for n=1:L
        if C_new(n) == 0
            continue;
        end
        m = m + 1;
        iPX(m) = (n-1)/L;
        iPY(m) = (C_new(n) - 1)*L/Nfft;
        iPW(m) = R_energy(C_new(n), n);
    end

    warning('off', 'polyfitweighted:RepeatedPointsOrRescale');
    p = polyfitweighted(iPX, iPY, degree_WPF, iPW);
    warning('on', 'polyfitweighted:RepeatedPointsOrRescale');
    PY = polyval(p, (0:L-1)/L);
%     pp = polyder(p);
%     PPY = polyval(pp, (0:L-1)/L);

    C_prev = C_new;
    KY_upper = round(PY*Nfft/L) + range_vec + 1;
    KY_upper = min(Nfft, max(1, KY_upper));
    KY_lower = round(PY*Nfft/L) - range_vec + 1;
    KY_lower = min(Nfft, max(1, KY_lower));
    for n=1:L
        k_upper = KY_upper(n);
        k_lower = KY_lower(n);
        [v, arg] = max(R_energy(k_lower:k_upper, n));
        k = arg + k_lower - 1;
        
        if v > 0
            C_new(n) = k;
        else
            C_new(n) = 0;
        end
    end
    
%     fprintf("ITER = %u\n", ITER);
%     CV = C_prev;
%     CV(CV == 0) = nan;
%     CV = CV';
%     CV2 = C_new;
%     CV2(CV2 == 0) = nan;
%     CV2 = CV2';
% 
%     figure;
%     imagesc(1:L, 1:Nfft, R_energy);
%     hold on;
%     plot(1:L, CV, 'r');
%     plot(1:L, CV2, 'g');
%     plot(1:L, KY_lower, 'm--');
%     plot(1:L, KY_upper, 'm--');
%     hold off;
%     set(gca,'ydir','normal');
%     colormap(flipud(gray));
%     axis square
%     colorbar;
%     title('PS');
%     pause;
end

% CV2 = C_new;
% CV2(CV2 == 0) = nan;
% CV2 = CV2';

% figure;
% imagesc(1:L, 1:Nfft, R_energy);
% hold on;
% plot(1:L, CV2, 'g');
% plot(1:L, KY_lower, 'm--');
% plot(1:L, KY_upper, 'm--');
% hold off;
% set(gca,'ydir','normal');
% colormap(flipud(gray));
% axis square
% colorbar;
% title('PS');
% pause;


%% curve completions and set rest of TF region to zero
R_out = R_energy;
for n=1:L
    k_upper = KY_upper(n);
    k_lower = KY_lower(n);
    R_out(k_lower:k_upper, n) = 0;
    if C_new(n) == 0
        [val, arg_max] = max(R_energy(k_lower:k_upper, n));
        if val > 0
            C_new(n) = k_lower + arg_max - 1;
        end
    end
end


end

