function [x] = FM_inverse(STFT, Fs, Nfft, g, cas)

%tfr  : STFT of signal x
%cas  : if 1, ...
%       if 2, ...
%       if 3, ...
%h    : filter
%x    : restored signal

[N_Y, L] = size(STFT);

% case without periodizing
x = zeros(L,1);
Lg = (length(g)-1)/2;
for n=1:L
    x(n) = Fs/g(Lg+1)*sum(STFT(:,n))/Nfft;
end
return;

% if (cas == 1)
% end

if (cas == 2)
    %cases with periodization
    Lg = (length(g)-1)/2;
    x = zeros(L,1);
    for n = 1:L
        ind  = n-Lg:n+Lg;
        if (n > Lg)&&(n <= L - Lg)
            x(n) = mean(STFT(:,ind).*exp(2*1i*pi*(0:N_Y-1)'*(n-ind)/N_Y)*...
            g(Lg+1+n-ind))/norm(g(Lg+1+n-ind))^2;
        else
            x(n) = mean((STFT(:,1+rem((ind-1)+L,L)).*exp(2*1i*pi*(0:N_Y-1)'*(n-ind)/N_Y))...
            *g(Lg+1+n-ind))/norm(g(Lg+1+n-ind))^2;
        end
    end
end

if (cas == 3)
    %case with periodization
    Lg = (length(g)-1)/2;
    x = zeros(L,1);
    for n = 1:L
        ind = n-Lg:n+Lg;
        if (n > Lg)&&(n <= L-Lg)
            x(n) = mean((STFT(:,ind).*exp(2*1i*pi*(0:N_Y-1)'*(n-ind)/N_Y))*...
            ones(length(Lg+1+n-ind),1))/sum(g(Lg+1+n-ind));
        else
            x(n) = mean((STFT(:,1+rem((ind-1)+L,L)).*exp(2*1i*pi*(0:N_Y-1)'*(n-ind)/N_Y))*...
            ones(length(Lg+1+n-ind),1))/sum(g(Lg+1+n-ind));
        end
    end
end
end

