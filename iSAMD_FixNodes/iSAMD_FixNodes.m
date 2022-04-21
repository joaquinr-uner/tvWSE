function [x, v_e,exitflag] = iSAMD_FixNodes(f,A,phi,D,vnv,vh,tn,mInterp,options)

if nargin<10
    options = optimoptions('lsqcurvefit','Algorithm','levenberg-marquardt','MaxFunctionEvaluations',400*2*12*9);
end

v0 = vh;
X = [A;phi];
N = length(f);
modelfun = @(v,X)(regresion_reform_fix(v,X,N,D,vnv,tn,mInterp));

warning('off')
%v_e = nlinfit([A;phi],f,modelfun,v0,options);
[v_e,~,~,exitflag] = lsqcurvefit(modelfun,v0,[A;phi],f,[],[],options);
warning('on')


x = modelfun(v_e,X);

end
function s = regresion_reform_fix(v,X,N,D,vnv,tn,mInterp)
I = length(D);
A = X(1,:);
phi = X(2,:);
s = A.*cos(2*pi*phi);

for i=1:D-1
    v_l = v(sum(vnv(1:i-1))+1+(i-1)*2:sum(vnv(1:i))+(i)*2);
    t_l = tn{i};
    A_l = v_l(1:vnv(i));
    %amp = ppval(pwch([1 t_l N],A_l,zeros(size(A_l))),1:length(s));
    amp = interp1(t_l,A_l,1:N,mInterp);
    e = v_l(end);
    gam = v_l(end-1);
    s = s + amp.*(cos(2*pi*e*phi)+ gam*sin(2*pi*e*phi));
    %v(2*(sum(vnv(1:i-1)))+1:2*sum(vnv(1:i))) = [t_l, A_l,gam,e];
end
%s = A.*s;
end
