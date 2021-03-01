function [modes, s_RM] = MR_simple(STFT, Fs, Nfft, g, Lg, K_lower, K_upper, Nr)

[N_Y, L] = size(STFT);
modes = zeros(Nr, L);

gamma_Vg = median(abs(real(STFT(:))))/0.6745;
STFT_TH = STFT.*(abs(STFT) > 2*gamma_Vg);

KYd = max(K_lower, 1);
KYu = min(K_upper, N_Y);
for n=1:L
    for p=1:Nr
        % if the size is of I_LC is too small, continue
        if KYd(p, n) == N_Y || KYu(p, n) == 1
            continue;
        end

        ensure_v = 1;
        % skip the case HIGH FREQ region below LOW FREQ
        if p == Nr || KYu(p + 1, n) <= KYd(p, n)
            ensure_v = 0;
        end

        % ensure disjoint intervals
        if ensure_v == 1 && KYu(p, n) >= KYd(p+1, n)
            mid = KYd(p+1, n) + round((KYu(p, n) - KYd(p+1, n))/2);
            KYu(p, n) = mid;
            KYd(p+1, n) = mid+1;
        end
    end
end

for p=1:Nr
    for n=1:L
        modes(p, n) = Fs/g(Lg+1)*sum(STFT_TH(KYd(p, n):KYu(p, n), n))/Nfft;
    end
end

s_RM = sum(modes, 1);

end

