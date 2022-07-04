function [a, b, Al, alp] = compute_ampl_coefs(ti_h,Ql,gamh,N,D,mInterp,nmn,B1)

if nargin<7
    mnm = 1;
end

B = cell(1,D-1);
alp = cell(1,D-1);
Al = zeros(D-1,N);
a = zeros(1,D);
b = zeros(1,D);
a(1) = sum(B1)/N;
for i=1:D-1
    ti = ti_h{i};
    if nmn == 0
        b1 = interp1(1:N,B1,ti,mInterp);
        alpi = Ql{i}./b1;
    else
        alpi = Ql{i};
    end
    Al(i,:) = interp1(ti,alpi,1:N,mInterp);
    a(i+1) = sum(Al(i,:),2)/N;
    Al(i,:) = Al(i,:)/a(i+1);
    b(i+1) = gamh(i).*a(i+1);
    alp{i} = alpi/a(i+1);
end
end

