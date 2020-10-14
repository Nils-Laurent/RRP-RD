function [rel_err_SNR_p] = relative_error_QM(s_vec, phip_vec, B_vec, p, snr_in)
rel_err_SNR_p = 0;

Nfft = 512;
L = length(s_vec(:, 1));
sigma_s = 1/sqrt(B_vec(p));

noise = randn(L,1)+1i*randn(L,1);
sn_LC = sigmerge(s_vec(:, p), noise, snr_in);
[g, Lg] = create_gaussian_window(L, Nfft, sigma_s);

std_lc_g = 1/(sqrt(2*pi)*sigma_s)*sqrt(1 + sigma_s^4*B_vec(p)^2);
std_lc_g = ceil(std_lc_g*Nfft/L);


[STFT, ~, ~, QM, Vxgp] = FM_operators(sn_LC, Nfft, g, Lg, sigma_s);

%% create LM matrices
% [STFT_LM] = LM_from_STFT(STFT);
% 
% gamma_Vg = median(abs(real(STFT(:))))/0.6745;
% gamma_Vxgp = median(abs(real(Vxgp(:))))/0.6745;
% 
% Filter_LM = STFT_LM.*(abs(STFT) > 3*gamma_Vg).*(abs(Vxgp) > 3*gamma_Vxgp);


% figure;
% imagesc(abs(STFT_LM.*(abs(STFT) > 3*gamma_Vg)));
% set(gca,'ydir','normal');
% axis square;
% colormap(flipud(gray));

%% test intersection
% TFR_th = STFT.*(abs(STFT) > 3*gamma_Vg);
% TFR_th2 = Vxgp.*(abs(Vxgp) > 3*gamma_Vxgp);
% TFR_prod = abs(TFR_th).*abs(TFR_th2);
% 
% % figure;
% imagesc(abs(TFR_th));
% set(gca,'ydir','normal');
% axis square;
% 
% figure;
% imagesc(abs(TFR_th2));
% set(gca,'ydir','normal');
% axis square;
% 
% figure;
% imagesc(TFR_prod);
% set(gca,'ydir','normal');
% axis square;

%% settings
SPEC = abs(STFT).^2;
IF_ref = round(phip_vec(:, p)*Nfft/L) + 1;

%% test zeros
th_1 = 10^(-6);

TFR_th = SPEC < th_1;

% TFR_fig = (SPEC + 3/2*mean(mean(SPEC))).*(1 - TFR_th);

% figure;
% imagesc(TFR_fig);
% set(gca,'ydir','normal');
% colormap(flipud(gray));
% axis square;

% return;

d_zeros = zeros(L, 1);
for n=1:L
    k0 = IF_ref(n);
    [I_gamma] = find(TFR_th(:, n));
    
    if isempty(I_gamma)
        d_zeros(n) = 1;
    else
        d_zeros(n) = exp(-min(abs(I_gamma - k0)));
    end
end

d_zeros = d_zeros/max(d_zeros);

%% relative error
TFR2 = zeros(size(STFT));

acc_n = [];
acc = [];

for n=1:L
    for k=2:(Nfft-1)
        if std_lc_g >= abs(k - IF_ref(n))
            if SPEC(k, n) > SPEC(k - 1, n)...
                && SPEC(k, n) > SPEC(k + 1, n)
                % local maximum at (k, n);
                
%                 acc = [acc, diff_CR^2/sum_CR^2];
                TFR2(k, n) = 1;
            end
        end
    end
    
    [v, argv] = max(SPEC(:, n).*TFR2(:, n));
    if v > 0
        k = argv;
        acc_n = [acc_n, n];

        diff_CR = real(QM(k, n)) - B_vec(p);
        sum_CR = abs(real(QM(k, n))) + B_vec(p);
        acc = [acc, diff_CR^2/sum_CR^2];
    end
end

%% display figures
[v_phi, phi] = max(TFR2.*SPEC);
diff_phi = max(0, abs(phi(2:end) - phi(1:end - 1)) - 1);

rq_vec = zeros(L, 1);
for n=1:L
    rq_vec(n) = real(QM(phi(n), n));
end
diff_rq = abs(rq_vec(2:end) - rq_vec(1:end-1));
diff_rq = diff_rq/max(diff_rq);

for n=1:L-1
    if v_phi(n) == 0
        diff_rq(n) = 1;
    end
end

if max(diff_phi > 0)
    diff_phi = diff_phi/max(diff_phi);
end

figure;
hold on;
plot(acc_n, acc);
plot(diff_rq, 'r--');
hold off;
title("diff rq");

% figure;
% hold on;
% plot(acc_n, acc);
% plot(diff_phi, 'r--');
% hold off;
% title("diff phi");
% 
% figure;
% hold on;
% plot(acc_n, acc);
% plot(d_zeros, "r--");
% hold off;
% title("dist zero");

figure;
imagesc(TFR2.*SPEC);
set(gca,'ydir','normal');
colormap(flipud(gray));
axis square;
title("max in std around IF");

% figure;
% imagesc(abs(Filter_LM));
% set(gca,'ydir','normal');
% axis square;
% colormap(flipud(gray));
% title("max > gamma");

rel_err_SNR_p = sqrt(mean(acc));

end

