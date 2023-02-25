function [x_r, v] = WSF(x,c,D,mode,b,F,sF)

global sigma
global fmax
if nargin<6
    switch mode
        case 1
            [F,sF] = STFT_Gauss(x,length(x),sigma,fmax);
        case 2
            [F1, sF] = STFT_Gauss(x,length(x),sigma,fmax);
            [F] = synchro_F(F1,0.5);
        case 3
            [F1, sF] = STFT_Gauss(x,length(x),sigma,fmax);
            F = threshold(F1);
        otherwise
    end
end
%c = ridge_ext(F,0.1,0.1,10,10);

y = zeros(size(x));

for i=1:length(y)
   binf = max([1,c(i)-b]);
   bup = min([size(F,1),c(i)+b]);
   y(i) = 1/max(sF)*sum(F(binf:bup,i)); 
end

A = abs(y);

phi = unwrap(angle(y))/(2*pi);


C = cosenos(A,phi,D);

%S = pinv(C'*C);

%S = C*S*C';

%nl = trace(S);

%tau1 = trace(S)/length(x);
%tau2 = trace(S*S)/length(x);

%g = (2*tau1-tau1^2/tau2)/(1-tau2)^2;
v = ((C'*C)\C')*x;

x_r = C*v;
end