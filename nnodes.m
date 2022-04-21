function vnv = nnodes(F,sF,c,del,D,nv,f)
    [~,N] = size(F);
    K = length(D);
    if nv == 0
        vnv = zeros(sum(D)-K,1);
        for k=1:K
            vnvk = zeros(D(k)-1,1);
            for i=2:D(k)
                [A,~] = extract_fundamentals(F,sF,i*c(k,:),del);
                Sh = abs(fft(A)).^2;
                Sh(1) = 0;
                acum = cumsum(Sh(1:N/2))/sum(Sh(1:N/2));
                indx = find(acum>0.8,1);
                vnvk(i-1) = ceil(2*f(indx))+1;
                %vnvk(i-1) = ceil(3*f(indx));
            end
            vnvk(vnvk<2) = 2;
            vnv(sum(D(1:k-1)-k+1)+1:sum(D(1:k))-k) = vnvk;
        end
        
    else
            vnv = nv*ones(sum(D)-K,1);
    end
end
