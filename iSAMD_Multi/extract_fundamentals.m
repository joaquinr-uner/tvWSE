function [A,phi] = extract_fundamentals(F,sF,c,b)

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
    end

A(k,:) = abs(y);
A(k,:) = sqrt(2)*A(k,:);
phi(k,:) = unwrap(angle(y))/(2*pi);
end
end