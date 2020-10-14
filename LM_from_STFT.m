function [STFT_LM] = LM_from_STFT(STFT)

ASTFT = abs(STFT);
STFT_LM = zeros(size(STFT));
[Nfft, L] = size(STFT);

for n=1:L
    for k=2:(Nfft-1)
        if ASTFT(k, n) > ASTFT(k - 1, n)...
            && ASTFT(k, n) > ASTFT(k + 1, n)
            STFT_LM(k, n) = STFT(k, n);
        end
    end
end


end

