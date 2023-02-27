function [U] = istct_fast(F,f,g,fpb)
% This function implements a fast computation of the inverse
% short-time cepstrum transform (iSTCT) for use in the
% calculation of the de-shape STFT
if nargin<4
    fpb = 0.03;
end

if nargin<3
    g = 0.3;
end

f = f+f(1);

[N, ~] = size(F);

if iscolumn(f)
    f = f';
end

E = exp(-2*pi*1i*repmat(f,[N 1]).*repmat(1./f',[1 N])); % Matriz tipo DFT, capaz con Vandermonde se puede hacer más rápido?

aux = abs(F).^g;
aux = aux - repmat(mean(aux),[N,1]);

D = designfilt('highpassiir','FilterOrder',1,'PassbandFrequency',fpb,'PassbandRipple',0.2,'SampleRate',1);
aux = filtfilt(D,aux);

U = E * aux;
end
