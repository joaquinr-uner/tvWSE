function [A,phi,y] = extract_harmonics(F,sF,c,b,h,D)
%%
% Harmonic decomposition of signal s with STFT F.
% Inputs:
%         F: STFT of signal s.
%         sF: sum of STFT filters. Used in the reconstruciton step.
%         c: fundamental ridge.
%         b: harmonic reconstruction half-width
%         h: ridge correction search width
%         D: number of harmonic compontents
% Outputs:
%         A: harmonic amplitudes of s
%         phi: harmonic phases of s

[M,N] = size(F);
A = zeros(D,N);
phi = zeros(D,N);

for k=1:D
    y = zeros(1,N);
    ck = ridge_correct(k*c,F,h);

    for i=1:length(y)
        binf = max([1,ck(i)-b]);
        bup = min([M,ck(i)+b]);
        y(i) = 1/max(sF)*sum(F(binf:bup,i)); 
        %figure(1)
        %plot(abs(F(binf:bup,i)))
        %pause(0.005)
    end

A(k,:) = abs(y);
phi(k,:) = phase(y)/(2*pi);
end
end
