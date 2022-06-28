function [x, v_e,exitflag] = iSAMD(f,A,phi,D,vnv,vh,mInterp,lb,ub,method,options)

if nargin<11 && method == 0
    options = statset;
    options.RobustWgtFun = 'cauchy';
    options.MaxIter = 3000;
else
    if nargin<11 && method == 1
        options = optimoptions('lsqcurvefit','Algorithm','levenberg-marquardt','MaxFunctionEvaluations',400*2*sum(vnv),'Display','off');
    end
end

v0 = vh;
X = [A;phi];
N = length(f);
modelfun = @(v,X)(regresion_reform(v,X,N,D,vnv,mInterp));

warning('off')
if method == 0
    fprintf('Using nlinfit...\n')
    v_e = nlinfit([A;phi],f,modelfun,v0,options);
else
    if method == 1
        fprintf('Using lsqcurvefit...\n')
        [v_e,~,~,exitflag] = lsqcurvefit(modelfun,v0,[A;phi],f,lb,ub,options);
    end
end
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
    %amp = interp1([1, t_l, N],A_l,1:N,mInterp);
    amp = interp1([0, t_l, 1],A_l,linspace(0,1,N),mInterp);
    e = v_l(end);
    gam = v_l(end-1);
    s = s + amp.*(cos(2*pi*e*phi)+ gam*sin(2*pi*e*phi));
    %subplot(2,1,i)
    %plot(amp)
end
end
