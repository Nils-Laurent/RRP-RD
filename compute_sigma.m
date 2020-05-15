%% signal definition
F = 1/4;
% F = 1;
L = 4096*F;
t = (0:L-1)'/L;

phi1 = 256*F*t+2700*F*(t.^2)/2;
phi2 = 768*F*t+3000*F*(t.^2)/2;
IF1 = 256*F+2700*F*t;
IF2 = 768*F+3000*F*t;
s1 = exp(2*1i*pi*phi1);
s2 = exp(2*1i*pi*phi2);
s_clean = s1 + s2;
% modes_LC = transpose([s1, s2]);
% IFs_LC = transpose([IF1, IF2]);
% sigma_LC = 0.0188;

Nfft = 512;

%% plot TFR
% sig_plot = 0.0188;
% [g, Lg] = create_gaussian_window(L, Nfft, sig_plot);
% [TFR, ~, ~, ~] = FM_operators(s_clean, Nfft, g, Lg, sig_plot);
% 
% figure;
% imagesc((0:L-1)/L, (0:Nfft-1)*L/Nfft, abs(TFR));
% set(gca,'ydir','normal');
% colormap(flipud(gray));
% axis square
% pause;

%% find sigma
sigma_set = 0.005:0.005:0.04;
SL = length(sigma_set);

[RE_vec] = sig_min(sigma_set, s_clean, L, Nfft);
[~, arg] = min(RE_vec);
fprintf('min at %u\n', arg);

% figure;
% plot(sigma_set, RE_vec);

% it 2
a = max(1, arg - 1);
b = min(SL, arg + 1);
sigma_set2 = sigma_set(a):0.001:sigma_set(b);

[RE_vec] = sig_min(sigma_set2, s_clean, L, Nfft);
[~, arg] = min(RE_vec);
fprintf('min at %u\n', arg);

% figure;
% plot(sigma_set2, RE_vec);

% it 3
a = max(1, arg - 1);
b = min(SL, arg + 1);
sigma_set3 = sigma_set2(a):0.0001:sigma_set2(b);

[RE_vec] = sig_min(sigma_set3, s_clean, L, Nfft);
[~, arg] = min(RE_vec);
fprintf('min at %u\n', arg);

% figure;
% plot(sigma_set3, RE_vec);

fprintf('sigma = %f\n', sigma_set3(arg));


