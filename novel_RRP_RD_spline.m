function [Cs, XCs, Qs, TFR_inter] = novel_RRP_RD_spline(TFR, QM, sigma_s, Nr, TOL)
scale = 8;
[Nfft, L] = size(TFR);
LNh = ceil(L/200);

% std_QM = 1/(sqrt(2*pi)*sigma_s)*sqrt( 1+sigma_s^4*real(QM(k, n))^2 );
sigma_I = 1/(sqrt(2*pi)*sigma_s);
range_I = ceil(3*sigma_I*Nfft/L);

%% compute matrix of LC based frequency sum
S_LM = novel_LM_sum(TFR, QM, sigma_s);
[~, S_LM_sorted] = sort(S_LM, 'descend');

%% compute stability information
fprintf("compute stable RRP\n");
[RRP_stable, RRP_E_stable] = novel_RRP_stable(TFR, S_LM, S_LM_sorted, QM, sigma_s, Nr, scale);

%% merge intersecting structures
fprintf("check if RRP are intersecting\n");
[RRP_inter, RRP_E_inter, ~] = novel_RRP_intersection(TFR, RRP_stable, RRP_E_stable, range_I, LNh);

%% spline

[~, id] = max(RRP_E_inter);
RRP_max = RRP_inter(id, :);
X = (1:L).*(RRP_max > 0);
X = X(X > 0);
XCs = X(1):X(end);
Y = RRP_max(RRP_max > 0);
[Qs, ~] = spaps((X-1)/L, Y, TOL);
Cs = round(fnval(Qs, (XCs - 1)/L));

% figure;
% imagesc(1:L, 1:Nfft, S_LM);
% set(gca,'ydir','normal');
% axis square
% colormap(flipud(gray));
% title("S LM");


TFR_stable = zeros(size(TFR));
for r=1:length(RRP_E_stable)
    for n=1:L
        kr = RRP_stable(r, n);
        if kr > 0
            TFR_stable(kr, n) = RRP_E_stable(r);
        end
    end
end

% figure;
% imagesc(1:L, 1:Nfft, TFR_stable);
% set(gca,'ydir','normal');
% axis square
% colormap(flipud(gray));
% title("TFR stable");

TFR_inter = zeros(size(TFR));
for r=1:length(RRP_E_inter)
    for n=1:L
        kr = RRP_inter(r, n);
        if kr > 0
            TFR_inter(kr, n) = RRP_E_inter(r);
        end
    end
end

% figure;
% imagesc(1:L, 1:Nfft, TFR_inter);
% set(gca,'ydir','normal');
% axis square
% colormap(flipud(gray));
% title("TFR inter");
% hold on;
% plot(a:b, Cs);
% hold off;

% figure;
% imagesc(1:L, 1:Nfft, abs(TFR));
% set(gca,'ydir','normal');
% axis square
% colormap(flipud(gray));
% title("TF regions");
% hold on;
% plot(1:L, KY_lower, 'b');
% plot(1:L, KY_upper, 'r');
% hold off;
% pause;

end

