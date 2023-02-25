function h1 = spectro_log(F,t,f)
%UNTITLED8 Summary of this function goes here
%   Detailed explanation goes here
if nargin<3
    imagesc(log10(abs(F).^2)), axis xy
else
    imagesc(t,f,log(abs(F).^2)), axis xy
end
colormap(1-gray)

end
