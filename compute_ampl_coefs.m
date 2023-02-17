function [a, b, Q, q] = compute_ampl_coefs(ti_h,Al,gamh,mInterp,nmn,B1)

if nargin<7
    mnm = 1;
end
D = length(gamh)+1;
N = length(B1);
B = cell(1,D-1);
q = cell(1,D-1);
Q = zeros(D-1,N);
a = zeros(1,D);
b = zeros(1,D);
a(1) = sum(B1)/N;
for i=1:D-1
    ti = ti_h{i};
    if nmn == 0
        b1 = interp1(1:N,B1,ti,mInterp);
        alpi = Al{i}./b1;
    else
        alpi = Al{i};
    end
    Q(i,:) = interp1(ti,alpi,1:N,mInterp);
    a(i+1) = sum(Q(i,:),2)/N;
    Q(i,:) = Q(i,:)/a(i+1);
    b(i+1) = gamh(i).*a(i+1);
    q{i} = alpi/a(i+1);
end
end

