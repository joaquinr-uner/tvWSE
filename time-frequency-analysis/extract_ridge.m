function A = extract_ridge(F,N,c)
% Extract the energy curve on the ridge c.
% Inputs:
%	  - F: Short-time Fourier Transform of signal x.
%	  - N: length of signal x.
%	  - c: extracted ridge. 
    [~,N] = size(F);
    A = zeros(1,N);
    for i=1:N
        A(i) = 1/N*abs(F(c(i),i)).^2;
    end
end