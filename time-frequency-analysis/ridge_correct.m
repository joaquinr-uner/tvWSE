function [cr] = ridge_correct (c,F,h,redun)
% correction of ridge location by band-limited maximum energy search
% Inputs:
%         c: harmonic ridge
%         F: STFT of signal s
%         h: ridge correction search width
%         redun: redundancy of the STFT
if nargin<4
   redun = 1; 
end
[K, N] = size(F);
cr = zeros(size(c));

for i=1:N
    binf = max([c(i)-h,1]);
    bup = min([c(i)+h,K]);
    I = binf:bup;
    [~,aux] = max(abs(F(I,i)).^2);
    if isempty(aux)
	cr(i) = c(i);
    else
	cr(i) = aux + binf - 1;
end

end