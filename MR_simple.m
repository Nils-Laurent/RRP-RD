function [modes, s_RM] = MR_simple(STFT, Nfft, XCs, g, Lg, K_lower, K_upper, Nr)

[N_STFT, L] = size(STFT);
modes = zeros(Nr, length(XCs));

KYd = K_lower;
KYu = K_upper;
for n=1:length(XCs)
    for p=1:Nr
        % if the size is of I_LC is too small, continue
        if KYd(p, n) == N_STFT || KYu(p, n) == 1
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
    for n=1:length(XCs)
        icol = XCs(n);
        modes(p, n) = L/g(Lg+1)*sum(STFT(KYd(p, n):KYu(p, n), icol))/Nfft;
    end
end

s_RM = sum(modes, 1);

end

