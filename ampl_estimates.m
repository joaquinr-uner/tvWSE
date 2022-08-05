function [a, Ah] = ampl_estimates(F,sF,c,del,D,vnv)
    [~,N] = size(F);
    K = length(D);
    Ah = zeros(sum(D)-K,N);
    for k=1:K
        for i=2:D(k)
            ck = ridge_correct(i*c(k,:),F,del,1);
            A = extract_ridge(F,N,ck);
            %[A,~] = extract_fundamentals(F,sF,ck,del);
            sh = A;
            idx_n = sum(D(1:k-1)-1)+i-1;
            nv = vnv(idx_n);
            inter = round(linspace(1,N,nv));
            idxs = sum(vnv(1:idx_n-1))+1:sum(vnv(1:idx_n));
            a(idxs)=sh(inter);
            Ah(i-1,:) = A;
        end
    end
end