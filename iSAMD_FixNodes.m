function [x, v_e] = iSAMD_FixNodes(f,A,phi,D,vnv,options,vh,tn,mInterp)
if isrow(f)
    f = f';
end

[I, N] = size(A);
if nargin<7
    options = statset;
    options.RobustWgtFun = 'cauchy';
    options.MaxIter = 100;
    mInterp = 'pchip';
end


if nargin<8
    C = zeros(2*sum(D),N);
    for i = 1:I
        C(2*sum(D(1:i-1))+1:2*sum(D(1:i)),:) = cosenos_marce(A(i,:),phi(i,:),D(i));
    end
    coef_est = (f*C')/(C*C');
    
    a = coef_est(2:D);
    b = coef_est(D+2:end);
    
    v0 = [];
    for i=1:length(vnv)
        inter = linspace(1,N,vnv(i));
        ax = inter(2:end-1);
        ayi = repmat(a(i),[1,vnv(i)]);
        bx = inter(2:end-1);
        byi = repmat(b(i),[1,vnv(i)]);
        v0 = [v0 ax, ayi, bx, byi];
        
    end
    v0 = [v0,2:D];
else
    v0 = vh;
end
X = [A;phi];

modelfun = @(v,X)(regresion_reform(v,N,X,D,vnv,tn,mInterp));

warning('off')
v_e = nlinfit([A;phi],f,modelfun,v0,options);
warning('on')


x = modelfun(v_e,X);
end
function s = regresion_reform(v,N,X,D,vnv,tn,mInterp)
        I = length(D);
        A = X(I,:);
        phi = X(2*I,:);
        s = A.*cos(2*pi*phi);
        for i=1:D-1
            v_l = v(sum(vnv(1:i-1))+1+(i-1)*2:sum(vnv(1:i))+(i)*2);
            t_l = tn{i};
            A_l = v_l(1:vnv(i));
            %amp = spline(t_l,A_l,t);
            amp = interp1(t_l,A_l,1:N,mInterp);
            e = v_l(end);
            gam = v_l(end-1);
            s = s + amp.*(cos(2*pi*e*phi)+ gam*sin(2*pi*e*phi));
        end 
        %s = A.*s;
        if isrow(s)
            s = s';
        end
end