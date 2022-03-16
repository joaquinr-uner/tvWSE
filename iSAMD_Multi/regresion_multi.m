function s = regresion_multi(v,t,X,D,vnv)
        K = length(D);
        s = 0;
        for k=1:K
            A = X(2*k-1,:);
            phi = X(2*k,:);
            s = s + A.*cos(2*pi*phi) + A.*sin(2*pi*phi);
            nvi = vnv(sum(D(1:k-1))-(k-2):sum(D(1:k)-1));
            vi = v(4*(sum(vnv(1:sum(D(1:k-1)-1)))-(sum(D(1:k-1)-1)))+1+sum(D(1:k-1))+(1-k):4*(sum(vnv(1:sum(D(1:k)-1)))-(sum(D(1:k)-1)))+sum(D(1:k-1))+(1-k));
            [alpha,beta] = compute_splines2(vi,t,nvi,D(k));
            e = v(4*(sum(vnv(1:sum(D(1:k)-1)))-(sum(D(1:k)-1)))+sum(D(1:k-1))+(1-k)+1:4*(sum(vnv(1:sum(D(1:k)-1)))-(sum(D(1:k)-1)))+sum(D(1:k-1))+(1-k)+D(k)-1);
            for j=1:D(k)-1
                s = s + A.*alpha(j,:).*cos(2*pi*e(j)*phi)...
                + A.*beta(j,:).*sin(2*pi*e(j)*phi);
            end
        end
end