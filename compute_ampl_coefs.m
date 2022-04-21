function [a, b, O, alp] = compute_ampl_coefs(ti_h,Ql,gamh,N,D,nv,B1,mInterp)

alp = zeros(D-1,nv);
O = zeros(D-1,N);
v = [];
a = zeros(1,D);
b = zeros(1,D);
a(1) = sum(B1)/N;
for i=1:D-1
    ti = ti_h{i};
    o = Ql{i};
    b1 = interp1(1:N,B1,ti,mInterp);
    O(i,:) = interp1(ti,o./b1,1:N,mInterp);
    a(i+1) = sum(O(i,:),2)/N;
    b(i+1) = gamh(i).*a(i+1);
    alp(i,:) = (o./b1)/a(i+1);
end
end

