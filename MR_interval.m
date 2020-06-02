function [K_lower, K_upper] = MR_interval(Cs, QM, Nfft, sigma_s)
    [Nr, L] = size(Cs);
    K_lower = zeros(size(Cs));
    K_upper = zeros(size(Cs));
    for p=1:Nr
        for n=1:L
            mod = real(QM(Cs(p, n), n));
            sigma_ILC = 1/(sqrt(2*pi)*sigma_s)*sqrt( 1+sigma_s^4*mod^2 );
            range = ceil(3*sigma_ILC*Nfft/L);
            K_lower(p, n) = Cs(p, n) - range;
            K_upper(p, n) = Cs(p, n) + range;
        end
        K_lower(p, :) = min(Nfft, max(1, K_lower(p, :)));
        K_upper(p, :) = min(Nfft, max(1, K_upper(p, :)));
    end
end

