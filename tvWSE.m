function [x_r, v_e,exitflag] = tvWSE(x,A,phi,D,vnv,vh,mInterp,lb,ub,method,ext,outn,options)
%% This function implements the adaptive time-varying wave-shape extraction (tvWSE) Algorithm
% detailed in Ruiz, J, et al. "Fully Adaptive Time-Varying Wave-Shape
% Model: Applications in Biomedical Signal Processing". 
% Inputs:
%         - x: Input Signal x.
%         - A: Amplitude Modulation (AM) of x.
%         - phi: phase modulation of x.
%         - D: number of harmonic components of the time-varying WSF.
%         - vnv: number of nodes per harmonic amplitude function.
%         - vh: initial coefficient vector.
%         - mInterp: Interpolation method for harmonic amplitude functions.
%         - lb: Lower bound for the coefficients in the lsqcurvefit algorithm.
%         - up: upper bound of the coefficients in the lsqcurvefit algorithm.
%         - method: Nonlinear regressiong algorithm. If == 0, statistical               
%                   non-linear curve fitting is performed. If ==1, least-square
%                   regression curve fitting is performed.
%         - ext: Indicates if the input signals has been extended.
%         - outn: Number of outer nodes in the harmonic amplitudes. Used
%                 when the signal is previously extended
%         - options: struct with optimization options for the curve fitting
%                    algorithms.
% Outputs:
%         - x_r: reconstructed signal obtained using the optimal
%         time-varying wave-shape coefficient vector.
%         - v_e: optimal time-varying wave-shape coefficient vector.
%         - exitflag: return the algorithm stop criterion for the curve
%                      fitting algorithm. See https://www.mathworks.com/help/optim/ug/lsqcurvefit.html
%                      for more details.

if nargin<11
    ext = 0;
end

if nargin<13 && method == 0
    options = statset;
    options.RobustWgtFun = 'cauchy';
    options.MaxIter = 3000;
else
    if nargin<13 && method == 1
        options = optimoptions('lsqcurvefit','Algorithm','levenberg-marquardt','MaxFunctionEvaluations',400*length(vh),'Display','off');
    end
end

v0 = vh;
X = [A;phi];
if ext == 1
    modelfun = @(v,X)(regresion_ext(v,X,D,vnv,mInterp,outn));
else
    modelfun = @(v,X)(regresion_reform(v,X,D,vnv,mInterp));
end
warning('off')
if method == 0
    fprintf('Using nlinfit...\n')
    v_e = nlinfit([A;phi],x,modelfun,v0,options);
else
    if method == 1
        fprintf('Using lsqcurvefit...\n')
        [v_e,~,~,exitflag] = lsqcurvefit(modelfun,v0,[A;phi],x,lb,ub,options);
    end
end
warning('on')


x_r = modelfun(v_e,X);
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
    amp = interp1([0, t_l, 1],A_l,linspace(0,1,N),mInterp);
    e = v_l(end);
    gam = v_l(end-1);
    s = s + amp.*(cos(2*pi*e*phi)+ gam*sin(2*pi*e*phi));
end
end

function s = regresion_ext(v,X,D,vnv,mInterp,outn)
I = length(D);
A = X(1,:);
phi = X(2,:);
s = A.*cos(2*pi*phi);
N = length(A);
for i=1:D-1
    v_l = v(2*(sum(vnv(1:i-1))+2*outn*(i-1))+1:2*(sum(vnv(1:i))+2*outn*i));
    t_l = v_l(1:vnv(i)+2*(outn-1));
    A_l = v_l(vnv(i)+2*(outn-1)+1:end-2);
    amp = interp1([0, t_l, 1],A_l,linspace(0,1,N),mInterp);
    e = v_l(end);
    gam = v_l(end-1);
    s = s + amp.*(cos(2*pi*e*phi)+ gam*sin(2*pi*e*phi));
end
end
