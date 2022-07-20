function [A,phi,y] = extract_fundamentals(F,sF,c,b)
% Función que extrae la amplitud y fase instantánea de la señal a
% partir de la reconstrucción en torno a la cresta fundamental en la STFT. 
[M,N] = size(F);
K = size(c,1);
A = zeros(K,N);
phi = zeros(K,N);
for k=1:K
    y = zeros(1,N);
    for i=1:length(y)
        binf = max([1,c(k,i)-b]);
        bup = min([M,c(k,i)+b]);
        y(i) = 1/max(sF)*sum(F(binf:bup,i)); 
        %figure(1)
        %plot(abs(F(binf:bup,i)))
        %pause(0.005)
    end

A(k,:) = abs(y);
phi(k,:) = phase(y)/(2*pi);
end
end