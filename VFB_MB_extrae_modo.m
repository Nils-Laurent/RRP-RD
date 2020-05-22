function [m,Fm] = VFB_MB_extrae_modo(F,c,q_,N,Nfft,sigma_opt)
     Fm = F;
     m = zeros(N,1);
     for  kkk = 1:N
        delta2 = ceil(3*(Nfft/N*(sqrt((1+sigma_opt^4*q_(c(kkk),kkk)^2)/(2*pi*sigma_opt^2)) )));
        m(kkk) = N/Nfft*sigma_opt*sum(F(max(1,c(kkk)-delta2):min(Nfft,c(kkk)+delta2),kkk));
        Fm(max(1,c(kkk)-delta2):min(Nfft,c(kkk)+delta2),kkk) = 0;
     end;
end