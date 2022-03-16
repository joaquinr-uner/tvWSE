function [x, v_e] = iSAMD_2step(f,t,A,phi,D,vnv,options,vh,ep)
    [I, N] = size(A);
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
    
    vha = [];
    for i=1:length(vnv)
        inter = linspace(t(1),t(end),vnv(i));
        ax = inter(2:end-1);
        ayi = repmat(a(i),[1,vnv(i)]);
        bx = inter(2:end-1);
        byi = repmat(b(i),[1,vnv(i)]);
        vha = [vha ax, ayi, bx, byi];
        
    end
    vhp = 2:D;
    %vh = [vha,vhp];
    
else
    vha = vh(1:4*(sum(vnv)-1)-D-1);
    vhp = vh(end-D+2:end);
end
X = [A;phi];

maxiter = options.MaxIter;
options.MaxIter = 1;

va = vha;
vp = vhp;

modelfun = @(v,X)(regresion(v,t,X,D,vnv));

modelfun1 = @(va,X)(regresion_over_amp(va,t,X,D,vnv));

modelfun2 = @(vp,X)(regresion_over_phi(vp,t,X,D,vnv));

x = zeros(size(f));
e = norm(x-f)/norm(f);
I = 0;

warning('off')
while e>ep && I<maxiter
    va = nlinfit({[A;phi],vp},f,modelfun1,va,options);
    
    vp = nlinfit({[A;phi],va},f,modelfun2,vp,options);

    v_e = [va, vp];
    x = modelfun(v_e,X);
    I = I+1;
    e = norm(x-f)/norm(f);
end

x = modelfun(v_e,X);
warning('on')
end

function s = regresion(v,t,X,D,vnv)
        I = length(D);
        A = X(I,:);
        phi = X(2*I,:);
        s = A.*cos(2*pi*phi) + A.*sin(2*pi*phi);
        [alpha,beta] = compute_splines2(v,t,vnv,D);
        for i=1:D-1
            e = v(4*(sum(vnv)-(D-1))+i);
            s = s + A.*alpha(i,:).*cos(2*pi*e*phi)...
            + A.*beta(i,:).*sin(2*pi*e*phi);
        end   
end

function s = regresion_over_amp(va,t,Xp,D,vnv)
        I = length(D);
        X = Xp{1};
        vp = Xp{2};
        A = X(I,:);
        phi = X(2*I,:);
        
        s = A.*cos(2*pi*phi) + A.*sin(2*pi*phi);
        [alpha,beta] = compute_splines2(va,t,vnv,D);
        for i=1:D-1
            e = vp(i);
            s = s + A.*alpha(i,:).*cos(2*pi*e*phi)...
            + A.*beta(i,:).*sin(2*pi*e*phi);
        end
end


function s = regresion_over_phi(vp,t,Xa,D,vnv)
        I = length(D);
        X = Xa{1};
        va = Xa{2};
        A = X(I,:);
        phi = X(2*I,:);
        
        s = A.*cos(2*pi*phi) + A.*sin(2*pi*phi);
        [alpha,beta] = compute_splines2(va,t,vnv,D);
        for i=1:D-1
            e = vp(i);
            s = s + A.*alpha(i,:).*cos(2*pi*e*phi)...
            + A.*beta(i,:).*sin(2*pi*e*phi);
        end
end