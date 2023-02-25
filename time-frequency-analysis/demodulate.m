function [s_demod,nonu_samp] = demodulate(t,s,iphi,ifrec,N)
    L = round(median(ifrec(0.1*N:0.9*N)));
    nonu_samp = iphi/L;
    nonu_samp = nonu_samp - nonu_samp(1);    
    esc = t(end)/nonu_samp(end);
    nonu_samp = nonu_samp*esc;
    s_demod = spline(nonu_samp,s,t);
end

