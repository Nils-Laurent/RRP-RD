function [g, Lh] = create_gaussian_window(Fs, Nfft, sigma_w)

prec = 10^(-3);
Lw =  Fs*sigma_w;
Lh = floor(Lw*sqrt(-log(prec)/pi))+1;
g = amgauss(2*Lh+1,Lh+1,Lw);

if 2*Lh + 1 > Nfft
%     sigma_w
    fprintf('[Warning] %f, 2*Lh+1 = %u > Nfft\n', sigma_w, 2*Lh + 1);
end

end