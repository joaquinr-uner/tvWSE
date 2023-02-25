function h1 = spectro(F,t,f)
% Compute and plot the Spectrogram from the Short-Time Fourier Transform.
% Inputs:
% 	  - F: Short-Time Fourier Transform.
%	  - t: time axis.
%	  - f: frequency axis.

if nargin<3
    imagesc(abs(F).^2), axis xy
else
    imagesc(t,f,abs(F).^2), axis xy
end
colormap(1-gray)

end

