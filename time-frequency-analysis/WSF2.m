function [x_r, x_r1, x_r2, v1, v2] = WSF2(x,c1,c2,D1,D2,mode,b,F,sF)

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

y1 = zeros(size(x));

for i=1:length(y1)
   binf = max([1,c1(i)-b]);
   bup = min([size(F,1),c1(i)+b]);
   y1(i) = 1/max(sF)*sum(F(binf:bup,i)); 
end

A1 = abs(y1);

phi1 = unwrap(angle(y1))/(2*pi);


y2 = zeros(size(x));

for i=1:length(y2)
   binf = max([1,c2(i)-b]);
   bup = min([size(F,1),c2(i)+b]);
   y2(i) = 1/max(sF)*sum(F(binf:bup,i)); 
end

A2 = abs(y2);

phi2 = unwrap(angle(y2))/(2*pi);

C1 = cosenos(A1,phi1,D1);

%v1 = ((C1'*C1)\C1')*x;

%x_r1 = C1*v1;
C2 = cosenos(A2,phi2,D2);

%v2 = ((C2'*C2)\C2')*x;

%x_r2 = C2*v2;

C = [C1 C2];

v = ((C'*C)\C')*x;

v1 = v(1:2*D1);

v2 = v(2*D1+1:end);

x_r1 = C1*v1;

x_r2 = C2*v2;

x_r = C*v;
end