function h1 = spectro(F,t,f)
%UNTITLED8 Summary of this function goes here
%   Detailed explanation goes here
if nargin<3
    imagesc(abs(F).^2), axis xy
else
    imagesc(t,f,abs(F).^2), axis xy
end
colormap(1-gray)

end

