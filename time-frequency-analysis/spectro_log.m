function h1 = spectro_log(F,t,f)
% Compute and plot the Log-spectrogram from the Short-Time Fourier Transform.
% Inputs:
% 	  - F: Short-Time Fourier Transform.
%	  - t: time axis.
%	  - f: frequency axis.
if nargin<3
    imagesc(log10(abs(F).^2)), axis xy
else
    imagesc(t,f,log(abs(F).^2)), axis xy
end
colormap(1-gray)

end
