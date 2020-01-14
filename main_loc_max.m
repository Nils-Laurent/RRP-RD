close all

%% signal definition
L = 4096;
t = (0:L-1)'/L;

A = 2048;
B = 3/4*L/2;
phi = A*t+B*cos(4*pi*t)/(4*pi);
s_clean = exp(2*1i*pi*phi);

%% ...
p2 = 4*pi*B;
sigma_Q = [0.0025, 1/((7/3)^(1/4)*sqrt(p2)), 1/sqrt(p2), 0.01];

J = 50;
Q = length(sigma_Q);

% for q=1:Q
%     [g, Lg] = create_gaussian_window(L, Nfft, sigma_Q(q));
%     TFR_clean = tfrstft(s_clean, Nfft, cas, g, Lg);
%     figure;
%     t_str = sprintf("\\sigma_%d = %f", q, sigma_Q(q));
%     imagesc((0:L-1)/L, (L/Nfft)*(1:Nfft), abs(TFR_clean));
%     set(gca,'ydir','normal');
%     title(t_str);
%     axis square
% end
% pause;

Mn = zeros(Q, J);
Mk = zeros(Q, J);

Nfft = 512;
Lg_X = Lg;
Lg_Y = round(Lg*Nfft/L);

for j=1:J
    if mod(j, 10) == 0
        fprintf("j = %d\n", j);
    end
    WGN = randn(L,1)+1i*randn(L,1);
    s_noise = sigmerge(s_clean, WGN, -10);
    for q=1:Q
        [g, Lg] = create_gaussian_window(L, Nfft, sigma_Q(q));
        [TFR_noise] = tfrstft(s_noise, Nfft, cas, g, Lg);
        g_lim = median(abs(real(TFR_noise(:))))/0.6745;
        Y = abs(TFR_noise);
        [Z, arg] = max(Y(:));
        g_lim = Z - g_lim;
        
        r = 1;
        while Z > g_lim
            [aY, aX] = ind2sub(size(Y),arg);
            Mk(q, j, r) = aY;
            Mn(q, j, r) = aX;
            
            x0 = max(1, aX-Lg_X);
            x1 = min(aX + Lg_X, L);
            y0 = max(1, aY-Lg_Y);
            y1 = min(aY + Lg_Y, Nfft);
            Y(y0:y1, x0:x1) = 0;
            [Z, arg] = max(Y(:));
            r = r + 1;
        end
    end
end

for q=1:Q
    Xs = Mn(q, :, :);
    Xs = Xs(:);
    Ys = Mk(q, :, :);
    Ys = Ys(:);
    t_str = sprintf("\\sigma_%d = %f", q, sigma_Q(q));
    [g, Lg] = create_gaussian_window(L, Nfft, sigma_Q(q));
    TFR_clean = tfrstft(s_clean, Nfft, cas, g, Lg);
    figure;
    imagesc((0:L-1)/L, (L/Nfft)*(1:Nfft), abs(TFR_clean));
    set(gca,'ydir','normal');
    title(t_str);
    axis square
    hold on;
    argX = floor(arg/Nfft);
    argY = mod(arg, Nfft);
    plot((Xs - 1)/L, (Ys - 1)*(L/Nfft), '+', 'MarkerSize', 10, 'MarkerEdgeColor','red');
    hold off;
end