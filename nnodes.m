function [vnv,alp] = nnodes(F,sF,c,del,D,nv,fs)
    [~,N] = size(F);
    K = length(D);
    if nv == 0
        vnv = zeros(sum(D)-K,1);
        alp = zeros(sum(D)-K,N);
        for k=1:K
            vnvk = zeros(D(k)-1,1);
            for i=2:D(k)
                [Al,~] = extract_fundamentals(F,sF,i*c(k,:),del);
                A = Al;
                %A = Al(0.1*N:0.9*N);
                Sh = abs(fft(A)).^2;
                Sh = Sh(1:floor(end/2));
                Sh(1) = 0;
                acum = cumsum(Sh)/sum(Sh);
                indx = find(acum>0.9,1);
                f = 0:fs/N:fs/2-fs/N;
                %f = 0:fs/(0.8*N):fs/2-fs/(0.8*N);
                vnvk(i-1) = 2*round(f(indx))+1;
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
