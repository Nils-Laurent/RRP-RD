function res = snr(s,n)
%SNR computes the SNR. s : signal and n: noise

rms = @(x) sqrt(mean(abs(x).^2));
res = 20 * log10(rms(s)/rms(n));

end

