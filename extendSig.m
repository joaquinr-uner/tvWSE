function [s_ext,spbk,spfw] = extendSig(s,phi,c,Np,method)
%%
% signal extension by ARIMA model forecasting
% Inputs: 
%         s: signal to extend
%         phi: fundamental phase function of s
%         c: number of cycles for ARIMA model fitting
%         Np: number of samples to forecast
%         method: choose between forward ('fw') or backward ('bw')
%         forecasting
% Outputs:
%         s_ext: extended signal
%         spbk: backward extended segment
%         spfw: forward extended segment

s = s(:);
phi = phi(:);

N = length(s);
wrapphi = wrapTo2Pi(2*pi*phi);
dwphi = diff(wrapphi);
plocs = find(dwphi<-5);
pw = diff(plocs);

%[~,plocs] = findpeaks(wrapTo2Pi(2*pi*phi));
%plocs = [1; plocs; N];

msplt = strsplit(method,'-');
spfw = [];
spbk = [];
if sum(ismember(msplt,'fw'))
    indfw = find(phi>phi(end)-c,1);
    Seas1 = round(median(pw(end-c+1:end)));
    estInd = (indfw:N)';
    if length(estInd)<Seas1+c %Check if enough samples are used for estimation
        estInd = N-Seas1-c-1:N;
    end

    Mdlfw = regARIMA('D',0,'Seasonality',Seas1,'MALags',c,'SMALags',Seas1,'Intercept',0);
    Mdlfwest = estimate(Mdlfw,s(estInd),'X',estInd,'Display','off');
    spfw = forecast(Mdlfwest,Np,'X0',(estInd),'Y0',s(estInd),'XF',(N+1:N+Np)');
    spfw = detrend(spfw);
end

if sum(ismember(msplt,'bw'))
    indbk = find(phi>phi(1)+c,1);
    Seas2 = round(median(pw(1:c)));
    sbk = flipud(s);
    estIndbk = (N-indbk+1:N)';
    if length(estIndbk)<Seas2+c %Check if enough samples are used for estimation
        estIndbk = N-Seas2-c-1:N;
    end

    Mdlbk = regARIMA('D',0,'Seasonality',Seas2,'MALags',c,'SMALags',Seas2,'Intercept',0);
    Mdlbkest = estimate(Mdlbk,sbk(estIndbk),'X',estIndbk,'Display','off');
    spbk = forecast(Mdlbkest,Np,'X0',(estIndbk),'Y0',sbk(estIndbk),'XF',(N+1:N+Np)');
    spbk = detrend(spbk);
end
    s_ext = [flipud(spbk); s; spfw];
    
end
