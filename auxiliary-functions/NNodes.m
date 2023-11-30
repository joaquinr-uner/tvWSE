function [vnv,alp] = NNodes(nv,D,F,sF,c,del,fs,thr)
%% Automatic node number estimation based on harmonic amplitude bandwidth.
% Inputs: 
% 	  - nv: Number of nodes. If nv == 0, adaptive number of nodes for each harmonic is computed.
%		Otherwise, number of nodes is set to nv for all harmonic.
%	  - D: Number of Harmonic Components.
%	  - F: Short-Time Fourier Transform of the signal x.
%	  - sF: Sum of filters of the STFT (used for vertical reconstruction).
%	  - c: fundamental ridge of the STFT.
%	  - fs: sampling frequency.
%	  - thr: threshold of the cumulative spectral energy.
% Outputs:
%	  - vnv: (D-1) length vector with the number of nodes for each harmonic.
%	  - alp: Rough estimates of the harmonic amplitudes obtained from F.

if nargin < 3
    K = length(D);
    vnv = nv*ones(sum(D)-K,1);
else
    if nargin < 8
        thr = 0.9;
    end
    [~,N] = size(F);
    K = length(D);
    if nv == 0
        vnv = zeros(sum(D)-K,1);
        alp = zeros(sum(D)-K,N);
        for k=1:K
            vnvk = zeros(D(k)-1,1);
            for i=2:D(k)
                ck = ridge_correct(i*c(k,:),F,del,1);
                Al = extract_ridge(F,N,ck);
                %[Al,~] = extract_fundamentals(F,sF,ck,del);
                A = Al(round(0.1*N)+1:round(0.9*N));
                Sh = abs(fft(A)).^2;
                Sh = Sh(1:floor(end/2));
                Sh(1) = 0;
                acum = cumsum(Sh)/sum(Sh);
                indx = find(acum>thr,1);
                Nc = length(A);
                f = 0:fs/Nc:fs/2-fs/Nc;
                %f = 0:fs/(0.8*N):fs/2-fs/(0.8*N);
                %vnvk(i-1) = 2*ceil(f(indx))+1;
                vnvk(i-1) = 2*indx+1;
                %vnvk(i-1) = ceil(3*f(indx));
                alp(sum(D(1:k-1)-k+1)+i-1,:) = Al;
            end
            vnvk(vnvk<2) = 2;
            vnv(sum(D(1:k-1)-k+1)+1:sum(D(1:k))-k) = vnvk;
        end


    else
        vnv = nv*ones(sum(D)-K,1);
    end
end
end
