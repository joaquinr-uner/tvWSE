function [C] = dict_iSAMD(gam,e,phi,D)

r = D-1;
%C = zeros(length(phi),2*r);
C = zeros(length(phi),r);
for i=1:r
%    C(:,i) = cos(2*pi*e(i)*phi);
%    C(:,r+i) = gam(i)*sin(2*pi*e(i)*phi);
     C(:,i) = cos(2*pi*e(i)*phi) + gam(i)*sin(2*pi*e(i)*phi);
end

end