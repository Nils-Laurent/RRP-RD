function [tfr] = tfrstft(x, N, cas, g, Lg)

[xrow,xcol] = size(x);
t = 1:xrow; %the time instant, we consider the time instant shitfed by a factor shift.
tfr = zeros (N,length(t));

if (cas == 1)
    %case without periodizing
    trans  = zeros(1,length(t));
    for icol=1:length(t)
        tau = -min([Lg,t(icol)-1]):min([Lg,xrow-t(icol)]);
        tfr(1:length(tau),icol) = x(t(icol)+tau,1).*g(Lg+1+tau);
        trans(icol)  = tau(1);
    end
    tfr=fft(tfr,N);
    A = exp(-2/N*pi*1i*(0:N-1)'*trans);
    tfr = tfr.*A;
end

if (cas == 2) || (cas == 3)
    %cases with periodization
    tau = -Lg:Lg;
    for icol = 1:length(t)
        if (t(icol) > Lg) && (t(icol) <= xrow - Lg)
            tfr(1:length(tau), icol) = x(t(icol) + tau,1).*g(Lg + 1 + tau);
        else
            tfr(1:length(tau), icol) = ...
                x(1+rem((t(icol)-1)+tau+xrow,xrow),1).*g(Lg+1+tau);
        end
    end
    tfr = fft(tfr, N);
    trans = Lg*ones(1, length(t));
    A = exp(2/N*pi*1i*(0:N-1)'*trans);
    tfr = tfr.*A;
end

end %tfrstft