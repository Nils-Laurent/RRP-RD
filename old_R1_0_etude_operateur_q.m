%% etude statistique de q

close all;
%% def

L = 4096;
t = (0:L-1)'/L;

A = L/32;
B_vec = [5*A, (5*A + L - 2*A)/2, (L - 2*A)];


s_vec = zeros(L, 3);
phip_vec = zeros(L, 3);
sigma_vec = zeros(3, 1);

for p=1:3
    phi = A*t + B_vec(p)*(t.^2)/2;
    s_vec(:, p) = exp(2*1i*pi*phi);
    sigma_vec(p) = 1/sqrt(B_vec(p));
    phip_vec(:, p) = A + B_vec(p)*t;
end

% figure;
% hold on;
% plot(t, phip_vec(:, 1));
% plot(t, phip_vec(:, 2));
% plot(t, phip_vec(:, 3));
% hold off;
% pause;

%% test

SNRs = -10;
L_snr = length(SNRs);
N_rep = 1;
RE_mean = zeros(3, L_snr);
RE_var = zeros(3, L_snr);
data_c = zeros(N_rep, 1);

for n_snr = 1:L_snr
    snr_in = SNRs(n_snr);
    for p=3:3
        fprintf("%u/%u, %u/3\n", n_snr, L_snr, p);
        for rep=1:N_rep
            data_c(rep) = relative_error_QM(s_vec, phip_vec, B_vec, p, snr_in);
        end
        RE_mean(p, n_snr) = mean(data_c);
        RE_var(p, n_snr) = std(data_c);
    end
end

% figure;
% hold on;
% plot(SNRs, RE_mean(1, :));
% plot(SNRs, RE_mean(2, :), '--');
% plot(SNRs, RE_mean(3, :), '.-');
% hold off;
% title("Relative error (mean)");
% 
% figure;
% hold on;
% plot(SNRs, RE_var(1, :));
% plot(SNRs, RE_var(2, :), '--');
% plot(SNRs, RE_var(3, :), '.-');
% hold off;
% title("Std");

