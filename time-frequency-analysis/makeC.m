function [C] = makeC(A,phi,D)

C = zeros(length(A),2*D);

for i=1:D
    C(:,i) = A.*cos(2*pi*i*phi);
    C(:,D+i) = A.*sin(2*pi*i*phi);
end

end