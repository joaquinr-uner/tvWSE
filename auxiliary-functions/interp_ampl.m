function A = interp_ampl(ti,alp,N,mInterp)
% Interpolation of optimal node locations and amplitudes.
% Inputs:
%	  - ti: cell structure with node locations for optimal nodes.
%	  - alp: cell structure with node amplitudes for optimal nodes.
%	  - N: length of signal.
%	  - mInterp: Interpolation method.
% Outputs:
%	  - A: Nx(D-1) matrix with harmonic amplitude functions.

if nargin<4
   mInterp = 'pchip';
end 

L = size(ti,2);
A = zeros(L,N);

    for i=1:L
        x_l = ti{i};
        y_l = alp{i};
        A(i,:) = interp1(x_l,y_l,1:N,mInterp);
    end
end