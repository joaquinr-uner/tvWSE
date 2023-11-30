function [a, b, Q, q] = compute_hafs(ti,alp,gamh,mInterp,nmn,B1)
%% Compute harmonic amplitude functions and coefficients.
% Perform normalization if necessary.
% Inputs:
%	  - ti: cell structure with node locations for optimal nodes.
%	  - alp: cell structure with node amplitudes for optimal nodes.
%	  - gamh: optimal phase-shift coefficients gamma_l.
%	  - mInterp: Interpolation method.
% 	  - nmn: Normalization flag. If equal to 1, perform normalization.
%	  - B1: Fundamental amplitude modulation. Required to perform normalization.
% Outputs:
%	  - a: Cosine amplitude coefficients.
%	  - b: Sine amplitude coefficients.
%	  - Q: Normalized harmonic amplitude functions.
%	  - q: cell structure with normalized node amplitudes.

if nargin<7
    mnm = 1;
end
D = length(gamh)+1;
N = length(B1);
B = cell(1,D-1);
q = cell(1,D-1);
Q = zeros(D-1,N);
a = zeros(1,D);
b = zeros(1,D);
a(1) = sum(B1)/N;
for i=1:D-1
    tii = ti{i};
    if nmn == 0
        b1 = interp1(1:N,B1,tii,mInterp);
        alpi = alp{i}./b1;
    else
        alpi = alp{i};
    end
    Q(i,:) = interp1(tii,alpi,1:N,mInterp);
    a(i+1) = sum(Q(i,:),2)/N;
    Q(i,:) = Q(i,:)/a(i+1);
    b(i+1) = gamh(i).*a(i+1);
    q{i} = alpi/a(i+1);
end
end

