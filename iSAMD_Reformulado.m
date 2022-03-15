function [x, v_e] = iSAMD_Reformulado(f,A,phi,D,vnv,t_corr,options,vh,mInterp)
if isrow(f)
    f = f';
end

[I, N] = size(A);
if nargin<6
    t_corr = 1;
end

if nargin<7
    options = statset;
    options.RobustWgtFun = 'cauchy';
    options.MaxIter = 100;
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

modelfun = @(v,X)(regresion_reform(v,N,X,D,vnv,t_corr,mInterp));

warning('off')
v_e = nlinfit([A;phi],f,modelfun,v0,options);
warning('on')


x = modelfun(v_e,X);
end
function s = regresion_reform(v,N,X,D,vnv,t_corr,mInterp)
I = length(D);
A = X(I,:);
phi = X(2*I,:);
s = A.*cos(2*pi*phi);
for i=1:D-1
    v_l = v(2*(sum(vnv(1:i-1)))+1:2*sum(vnv(1:i)));
    t_l = v_l(1:vnv(i)-2);
    if t_corr
        if t_l(1) <= 1
            %t_l(1) = round((t_l(2)+1)/2);
            t_l(1) = N/vnv(i);
        end
        
        if t_l(end) >= N
            %t_l(end) = round((N+t_l(end-1))/2);
            t_l(end) = N - N/vnv(i);
        end
        for j=2:length(t_l)-1
            if t_l(j)<=1 || t_l(j) >=N || t_l(j) == t_l(j+1)
                %t_l(j) = round((t_l(j-1)+t_l(j+1))/2);
                t_l(j) = j*N/vnv(i);
            end
        end
    end
    t_l = sort(t_l);
    A_l = v_l(vnv(i)-1:2*vnv(i)-2);
    %amp = ppval(pwch([1 t_l N],A_l,zeros(size(A_l))),1:length(s));
    amp = interp1([1, t_l, N],A_l,1:N,mInterp);
    e = v_l(end);
    gam = v_l(end-1);
    s = s + amp.*(cos(2*pi*e*phi)+ gam*sin(2*pi*e*phi));
    %v(2*(sum(vnv(1:i-1)))+1:2*sum(vnv(1:i))) = [t_l, A_l,gam,e];
end
%s = A.*s;
if isrow(s)
    s = s';
end
end