function [x] = itfrstft(tfr,cas,g)

%tfr  : STFT of signal x
%cas  : if 1, ...
%       if 2, ...
%       if 3, ...
%h    : filter
%x    : restored signal

[N,xrow] = size(tfr);
if (cas == 1)
    %case without periodizing
    x = zeros(xrow,1);
    Lg = (length(g)-1)/2;
    for icol=1:xrow
        x(icol) = 1/g(Lg+1)*mean(tfr(:,icol));
    end
end

if (cas == 2)
    %cases with periodization
    Lg = (length(g)-1)/2;
    x = zeros(xrow,1);
    for i = 1:xrow
        ind  = i-Lg:i+Lg;
        if (i > Lg)&&(i <= xrow - Lg)
            x(i) = mean(tfr(:,ind).*exp(2*1i*pi*(0:N-1)'*(i-ind)/N)*...
            g(Lg+1+i-ind))/norm(g(Lg+1+i-ind))^2;
        else
            x(i) = mean((tfr(:,1+rem((ind-1)+xrow,xrow)).*exp(2*1i*pi*(0:N-1)'*(i-ind)/N))...
            *g(Lg+1+i-ind))/norm(g(Lg+1+i-ind))^2;
        end
    end
end

if (cas == 3)
    %case with periodization
    Lg = (length(g)-1)/2;
    x = zeros(xrow,1);
    for i = 1:xrow
        ind = i-Lg:i+Lg;
        if (i > Lg)&&(i <= xrow-Lg)
            x(i) = mean((tfr(:,ind).*exp(2*1i*pi*(0:N-1)'*(i-ind)/N))*...
            ones(length(Lg+1+i-ind),1))/sum(g(Lg+1+i-ind));
        else
            x(i) = mean((tfr(:,1+rem((ind-1)+xrow,xrow)).*exp(2*1i*pi*(0:N-1)'*(i-ind)/N))*...
            ones(length(Lg+1+i-ind),1))/sum(g(Lg+1+i-ind));
        end
    end
end

end