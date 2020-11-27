function [STFT, omega, omega2, q, STFT_xgp, tau] = FM_operators(s, Fs, Nfft, g, Lg, sigma_s)
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

STFT_xgp   = zeros(Nfft,L);
% STFT_x2g   = zeros(Nfft,L);

ft = (0:Nfft-1)';

% Window related variables
t0  = ((0:2*Lg)'-Lg)/Fs;
gp = (-2*pi/sigma_s^2*t0).*g;
gpp = (-2*pi/sigma_s^2 + 4*pi^2/sigma_s^4*(t0.^2)).*g;

for n=1:L
    % STFT, window g
    time_inst = -min([Lg,n-1]):min([Lg,L-n]);
    vg = fft(s(n+time_inst).*g(Lg+time_inst+1),Nfft)/Fs;

    % STFT, window xg
    vxg = fft(s(n+time_inst).*(time_inst)'/Fs.*g(Lg+time_inst+1),Nfft)/Fs;
    
%     STFT_x2g(:, n) = fft(s(n+time_inst).*(time_inst.^2)'/L.*g(Lg+time_inst+1),Nfft)/L;

    % operator Lx (dtau)
    tau(:,n)  = vxg./vg;

    % STFT, window gp
    vgp = fft(s(n+time_inst).*gp(Lg+time_inst+1),Nfft)/Fs;


    % operator omega
    omega(:,n) = Fs/Nfft*ft-real(vgp/2/1i/pi./vg);

    % STFT, window gppsip(1:Lg-1) = diff(psi(1:Lg))*L^2/Nfft;p
    vgpp  = fft(s(n+time_inst).*gpp(Lg+time_inst+1),Nfft)/Fs;

    % STFT, windox xgp
    vxgp  = fft(s(n+time_inst).*(time_inst)'/Fs.*gp(Lg+time_inst+1),Nfft)/Fs;
    STFT_xgp(:, n) = vxgp;

    % computation of the two different omega

    q(:,n) = 1/2/1i/pi*(vgpp.*vg-vgp.^2)./(vxg.*vgp-vxgp.*vg);

    % new omega2
    omega2(:,n) = omega(:,n) - real(q(:,n)).*real(tau(:,n))...
                              + imag(q(:,n)).*imag(tau(:,n));

	% Storing STFT
    STFT(:,n) = vg.*(exp(2*1i*pi*min(Lg,n-1)*ft/Nfft));
end

end