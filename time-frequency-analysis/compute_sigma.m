function [sigma, fh, b] = compute_sigma(x,fhmode)
%% 
% Window parameter sigma automatic estimation from the power spectrum of
% signal x
% Inputs: 
%         x: signal under analysis
%         fhmode: average fundamental frequency estimation method. chose
%         from
%               '1': fh as Power Spectrum peak
%               '2': fh as the inverse of the fundamental period from the
%               cesptrum.
% Outputs:
%         sigma: window parameter sigma
%         fh: estimated signal dominant frequency
%         b: window half-support in frequency domain

if nargin<2
    fhmode = 1;
end

N = length(x);

X = fft(x);
X = X(1:floor(N/2));

switch fhmode
    case 1
        % fh as peak of the power spectrum
        [~,indf] = max(abs(X).^2);

        fh = indf+1;
    case 2
        % fh as inverse of fundamental period from complex cepstrum
        C = cceps(x);
        [~,idth] = max(C(floor(0.2*N):length(X)));
        th = idth+floor(0.2*N);
        fh = round(1/th*N);
    otherwise
        fprintf('Error during fh estimation')
end

sigma = log(100)/(5*1/fh*N)^2;

b = round(3/pi*sqrt(sigma/2)*N);

end