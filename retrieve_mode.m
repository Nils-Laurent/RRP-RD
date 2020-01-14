function [phipE_o1, phipE_o2, phippE] = retrieve_mode(s, Nfft, g, Lg, sigma_s, ridge)
%% retrieve_phase : retrieve signal phase with its' first and second derivatives
%
% INPUTS:
% s: real or complex signal
%    has to be a single mode
%    the phase of the signal has to be 0 at time t=0
% Nfft: frequency bins
% g: gaussian window s.t. g(x) = exp(-pi*x^2/sigma_s^2)
% Lg: gaussian window length
% sigma_s: gaussian coefficient
%
% OUTPUTS:
% STFT: STFT of the signal
% phi: phase
% omega2: phase second order first derivative
% phipp: phase second derivative

s = s(:);

L = length(s);
omega  = zeros(Nfft,L);
omega2 = zeros(Nfft,L);
q      = zeros(Nfft,L);
tau    = zeros(Nfft,L);
STFT   = zeros(Nfft,L);

ft = (0:Nfft-1)';

% Window related variables
t0  = ((0:2*Lg)'-Lg)/L;
gp = (-2*pi/sigma_s^2*t0).*g;
gpp = (-2*pi/sigma_s^2 + 4*pi^2/sigma_s^4*(t0.^2)).*g;

for n=1:L
    % STFT, window g
    time_inst = -min([Lg,n-1]):min([Lg,L-n]);
    vg = fft(s(n+time_inst).*g(Lg+time_inst+1),Nfft)/L;

    % STFT, window xg
    vxg = fft(s(n+time_inst).*(time_inst)'/L.*g(Lg+time_inst+1),Nfft)/L;

    % operator Lx (dtau)
    tau(:,n)  = vxg./vg;

    % STFT, window gp
    vgp = fft(s(n+time_inst).*gp(Lg+time_inst+1),Nfft)/L;


    % operator omega
    omega(:,n) = L/Nfft*ft-real(vgp/2/1i/pi./vg);

    % STFT, window gppsip(1:Lg-1) = diff(psi(1:Lg))*L^2/Nfft;p
    vgpp  = fft(s(n+time_inst).*gpp(Lg+time_inst+1),Nfft)/L;

    % STFT, windox xgp
    vxgp  = fft(s(n+time_inst).*(time_inst)'/L.*gp(Lg+time_inst+1),Nfft)/L;

    % computation of the two different omega

    q(:,n) = 1/2/1i/pi*(vgpp.*vg-vgp.^2)./(vxg.*vgp-vxgp.*vg);

    % new omega2
    omega2(:,n) = omega(:,n) - real(q(:,n)).*real(tau(:,n))...
                              + imag(q(:,n)).*imag(tau(:,n));

	% Storing STFT
    STFT(:,n) = vg.*(exp(2*1i*pi*min(Lg,n-1)*ft/Nfft));
end

%% compute estimates
phipE_o1 = zeros(L, 1);
phipE_o2 = zeros(L, 1);
phippE = zeros(L, 1);

for n = 1:L
    phipE_o1(n) = omega(ridge(n), n);
    phipE_o2(n) = omega2(ridge(n), n);
    phippE(n) = real(q(ridge(n), n));
end

% for n = 1:L
%     k = ridge(n);
%     th = 1/sqrt(2*pi)*sqrt(1/sigma_s^2 + sigma_s^2*real(phippE(n))^2);
%     th = round(th*Nfft/L);
%     lower = max(1, k - th);
%     upper = min(Nfft, k + th);
%     phippE_median(n) = median(real(q(lower:upper, n)));
% %     AmpE(n) = abs(STFT(ridge(n), n))/sigma_s*(1 + phippE(n)^2*sigma_s^4)^(1/4);
% 
% %     if (n == floor(L/2))
% %         V_test = zeros(size(ft));
% %         V_test(lower:upper) = phippE(n);
% %     end
%  end

% figure;
% hold on;
% plot((0:Nfft-1)*L/Nfft, real(q(:, floor(L/2))));
% plot((0:Nfft-1)*L/Nfft, V_test);
% hold off;
% pause;

%% extrapolation of phipE
%phipex = 200+1000*(0:L-1)/L;
% 
% phipE(1:Lg) = pchip(Lg:3*Lg,phipE(Lg:3*Lg),1:Lg);
% phipE(L-Lg:L) = pchip(L-2*Lg:L-Lg,phipE(L-2*Lg:L-Lg),L-Lg:L);

%% finite differences of phi'
% phippE(1:Lg-1) = diff(phipE(1:Lg))*L;
% phippE(L-Lg+1:L) = diff(phipE(L-Lg:L))*L;

end