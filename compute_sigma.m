close all;

%% signal definition
load fig1-observed-Livingston.txt
t = fig1_observed_Livingston(:,1);
t = t(:);
s = fig1_observed_Livingston(:,2);
s = s(:);
s = hilbert(s);
L = length(s);
% sigma_s = 0.05;
Nfft = 1024;
s_clean = s;

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
sigma_set = 0.005:0.005:0.1;
SL = length(sigma_set);

[RE_vec] = sig_min(sigma_set, s_clean, L, Nfft);
[~, arg] = min(RE_vec);
fprintf('min at %u\n', arg);

% figure;
% plot(sigma_set, RE_vec);
% pause;

% it 2
a = max(1, arg - 1);
b = min(length(sigma_set), arg + 1);
sigma_set2 = sigma_set(a):0.001:sigma_set(b);

[RE_vec] = sig_min(sigma_set2, s_clean, L, Nfft);
[~, arg] = min(RE_vec);
fprintf('min at %u\n', arg);

% figure;
% plot(sigma_set2, RE_vec);
% pause;

% it 3
a = max(1, arg - 1);
b = min(length(sigma_set2), arg + 1);
sigma_set3 = sigma_set2(a):0.0001:sigma_set2(b);

[RE_vec] = sig_min(sigma_set3, s_clean, L, Nfft);
[~, arg] = min(RE_vec);
fprintf('min at %u\n', arg);

% figure;
% plot(sigma_set3, RE_vec);
% pause;

fprintf('sigma = %f\n', sigma_set3(arg));


