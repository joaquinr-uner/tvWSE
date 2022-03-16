function [x, v_e, modes] = iSAMD_multi(f,t,A,phi,D,vnv,options,vh)
    [I, N] = size(A);
if nargin<7
    options = statset;
    options.RobustWgtFun = 'cauchy';
    options.MaxIter = 100;
end

if nargin<8

    C = construct_dct(A,phi,D);
    coef_est = (f*C')/(C*C');
    
    v0 = [];
    for i = 1:I
        g = coef_est(2*sum(D(1:i-1))+1:2*sum(D(1:i)));
        a = g(2:D(i));
        b = g(D(i)+2:end);
        vnv = Vnv(2*sum(D(1:i-1)+1:2*sun(D(1:i))));
        v = zeros(1,4*(sum(vnv)-(D(i)-1))+D(i)-1);
        for j=1:length(vnv)
            inter = linspace(t(1),t(end),vnv(j));
            ax = inter(2:end-1);
            ayj = repmat(a(j),[1,vnv(j)]);
            bx = inter(2:end-1);
            byj = repmat(b(j),[1,vnv(j)]);
            aux = [ax, ayj, bx, byj];
            v0(4*(sum(vnv(1:j-1))-(j-1))+1:4*(sum(vnv(1:j))-j)) = aux;
        end
        v0(4*(sum(vnv)-(D-1))+1:end) = 2:D;
    end
else
    v0 = vh;
end
X = [A;phi];

modelfun = @(v,X)(regresion_multi(v,t,X,D,vnv));

warning('off')
v_e = nlinfit([A;phi],f,modelfun,v0,options);
warning('on')


x = modelfun(v_e,X);

modes = zeros(length(D),N);
for k=1:length(D)
    Di = D(k);
    nvi = vnv(sum(D(1:k-1))-(k-2):sum(D(1:k)-1));
    vi = v_e(4*(sum(vnv(1:sum(D(1:k-1)-1)))-(sum(D(1:k-1)-1)))+1+sum(D(1:k-1))+(1-k):4*(sum(vnv(1:sum(D(1:k)-1)))-(sum(D(1:k)-1)))+sum(D(1:k-1))+(1-k));
    e = v_e(4*(sum(vnv(1:sum(D(1:k)-1)))-(sum(D(1:k)-1)))+sum(D(1:k-1))+(1-k)+1:4*(sum(vnv(1:sum(D(1:k)-1)))-(sum(D(1:k)-1)))+sum(D(1:k-1))+(1-k)+D(k)-1);
    modes(k,:) = regresion_multi([vi,e],t,[A(k,:);phi(k,:)],Di,nvi);
end
end

function s = regresion_multi(v,t,X,D,vnv)
        K = length(D);
        s = 0;
        for k=1:K
            A = X(k,:);
            phi = X(K+k,:);
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