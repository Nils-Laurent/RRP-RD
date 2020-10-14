%% etude de la propotion des maximum autour de l'IF

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

%% test

SNRs = -10:2:10;
L_snr = length(SNRs);
N_rep = 50;
prop_mean = zeros(3, L_snr);
prop_var = zeros(3, L_snr);
data_c = zeros(N_rep, 1);

for n_snr = 1:L_snr
    snr_in = SNRs(n_snr);
    for p=1:3
        fprintf("%u/%u, %u/3\n", n_snr, L_snr, p);
        for rep=1:N_rep
            [c_IF, c_ext] = count_IF_and_ext(s_vec, phip_vec, B_vec, p, snr_in);
            data_c(rep) = c_IF/(c_IF + c_ext);
        end
        prop_mean(p, n_snr) = mean(data_c);
        prop_var(p, n_snr) = std(data_c);
    end
end

figure;
hold on;
plot(SNRs, prop_mean(1, :));
plot(SNRs, prop_mean(2, :), '--');
plot(SNRs, prop_mean(3, :), '.-');
hold off;
title("Relative error (mean)");

figure;
hold on;
plot(SNRs, prop_var(1, :));
plot(SNRs, prop_var(2, :), '--');
plot(SNRs, prop_var(3, :), '.-');
hold off;
title("Std");
