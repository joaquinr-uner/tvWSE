function [a, Ah] = ampl_estimates(F,sF,c,del,D,vnv,N,Next)
K = length(D);
Ah = zeros(sum(D)-K,Next);
for k=1:K
    for i=2:D(k)
        [A,~] = extract_fundamentals(F,sF,i*c(k,:),del);
        sh = A;
        idx_n = sum(D(1:k-1)-1)+i-1;
        nv = vnv(idx_n);
        ti = ((Next-N+1)/2)/Next;
        te = 1-((Next-N)/2)/Next;
        inter = zeros(1,nv+2);
        inter(2:end-1) = round(linspace(ceil(ti*Next),floor(te*Next),nv));
        inter(1) = 1;
        inter(end) = N;
        idxs = sum(vnv(1:idx_n-1)+2)+1:sum(vnv(1:idx_n)+2);
        a(idxs)=sh(inter);
        Ah(i-1,:) = A;
    end
end
end