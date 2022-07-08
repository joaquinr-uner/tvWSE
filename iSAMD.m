function [x, v_e,exitflag] = iSAMD(f,A,phi,D,vnv,vh,mInterp,lb,ub,method,ext,options)
if nargin<11
    ext = 0;
end

if nargin<12 && method == 0
    options = statset;
    options.RobustWgtFun = 'cauchy';
    options.MaxIter = 3000;
else
    if nargin<12 && method == 1
        options = optimoptions('lsqcurvefit','Algorithm','levenberg-marquardt','MaxFunctionEvaluations',400*length(vh),'Display','off');
    end
end

v0 = vh;
X = [A;phi];
if ext == 1
    modelfun = @(v,X)(regresion_ext(v,X,D,vnv,mInterp));
else
    modelfun = @(v,X)(regresion_reform(v,X,D,vnv,mInterp));
end
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
function s = regresion_reform(v,X,D,vnv,mInterp)
I = length(D);
A = X(1,:);
phi = X(2,:);
s = A.*cos(2*pi*phi);
N = length(A);
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

function s = regresion_ext(v,X,D,vnv,mInterp)
I = length(D);
A = X(1,:);
phi = X(2,:);
s = A.*cos(2*pi*phi);
N = length(A);
for i=1:D-1
    v_l = v(2*(sum(vnv(1:i-1))+2*(i-1))+1:2*(sum(vnv(1:i))+2*i));
    t_l = v_l(1:vnv(i));
    A_l = v_l(vnv(i)+1:2*vnv(i)+2);
    %amp = interp1([1, t_l, N],A_l,1:N,mInterp);
    amp = interp1([0, t_l, 1],A_l,linspace(0,1,N),mInterp);
    e = v_l(end);
    gam = v_l(end-1);
    s = s + amp.*(cos(2*pi*e*phi)+ gam*sin(2*pi*e*phi));
    %subplot(2,1,i)
    %plot(amp)
    %title(num2str(e))
    %pause(0.1)
end
end
