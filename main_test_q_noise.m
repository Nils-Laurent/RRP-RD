close all

%% signal definition
L = 4096;
t = (0:L-1)'/L;

A = 45;
B = 4000;
phi = A*t+B*(t.^2)/2;
s_clean = exp(2*1i*pi*phi);
phip = A + B*t;
phipp = B*ones(L, 1);

%% max values detection
sigma = 1/sqrt(B);

Nfft = 512;
cas = 1;

WGN = randn(L,1)+1i*randn(L,1);
s_noise = sigmerge(s_clean, WGN, -10);

[g, Lg] = create_gaussian_window(L, Nfft, sigma);
%[TFR_noise] = tfrstft(s_noise, Nfft, cas, g, Lg);

%% 2nd order computation
[TFR_noise, omega, omega2, q] = FM_operators(s_noise, Nfft, g, Lg, sigma);
q_p2 = zeros(size(phipp));

gamma = median(abs(real(TFR_noise(:))))/0.6745;
xGamma = [];
yGamma = [];
max_TFR = max(abs(TFR_noise(:)));
TH_Gamma = 2.14*gamma;
% if (max_TFR - 3*gamma) > 3*gamma
%     fprintf("TH = max -3gamma\n")
%     TH_Gamma = max(3*gamma, max_TFR - 3*gamma);
% end

m = 1;
for n=1:L
    k = round(phip(n)*Nfft/L)+1;
    q_p2(n) = real(q(k, n));
    if (abs(TFR_noise(k, n)) > TH_Gamma)
        xGamma(m) = (n-1)/L;
        yGamma(m) = q_p2(n);
        m = m+1;
    end
end

% figure;
% imagesc((0:L-1)/L, (L/Nfft)*(1:Nfft), abs(TFR_noise));
% set(gca,'ydir','normal');
% title("noise");
% axis square

figure;
hold on;
%plot((0:L-1)/L, q_p2);
%plot((0:L-1)/L, phipp);
%plot((0:L-1)/L, zeros(size(phipp)));
plot((0:L-1)/L, (phipp - q_p2)/B);
plot((0:L-1)/L, zeros(size(phipp)));
plot(xGamma, (B - yGamma)/B, 'o');
hold off;