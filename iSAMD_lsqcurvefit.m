function [x, v_e,exitflag] = iSAMD_lsqcurvefit(f,A,phi,D,vnv,vh,mInterp,lb,ub,options)

if nargin<10
    options = optimoptions('lsqcurvefit','Algorithm','levenberg-marquardt','MaxFunctionEvaluations',400*2*12*9);
end

v0 = vh;
X = [A;phi];
N = length(f);
modelfun = @(v,X)(regresion_reform(v,X,N,D,vnv,mInterp));

warning('off')
%v_e = nlinfit([A;phi],f,modelfun,v0,options);
[v_e,~,~,exitflag] = lsqcurvefit(modelfun,v0,[A;phi],f,lb,ub,options);
warning('on')


x = modelfun(v_e,X);

end
function s = regresion_reform(v,X,N,D,vnv,mInterp)
I = length(D);
A = X(1,:);
phi = X(2,:);
s = A.*cos(2*pi*phi);

for i=1:D-1
    v_l = v(2*(sum(vnv(1:i-1)))+1:2*sum(vnv(1:i)));
    t_l = v_l(1:vnv(i)-2);
    A_l = v_l(vnv(i)-1:2*vnv(i)-2);
    %amp = ppval(pwch([1 t_l N],A_l,zeros(size(A_l))),1:length(s));
    amp = interp1([1, t_l, N],A_l,1:N,mInterp);
    e = v_l(end);
    gam = v_l(end-1);
    s = s + amp.*(cos(2*pi*e*phi)+ gam*sin(2*pi*e*phi));
    %v(2*(sum(vnv(1:i-1)))+1:2*sum(vnv(1:i))) = [t_l, A_l,gam,e];
end
%s = A.*s;
end
