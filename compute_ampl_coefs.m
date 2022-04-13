function [a, b, alp, v] = compute_ampl_coefs(ti_h,Al,gamh,eh,N,D,vlr,B1)

alp = zeros(D-1,N);
a = sum(Al,2)/N;
b = gamh.*a';
v = [];
for i=1:D-1
    alp(i,:) = (Al(i,:)./B1)/a(i);
    aux = [ti_h{i} alp(i,ti_h{i}) a(i) b(i) eh(i)];
    v = [v aux];
end
a1 = sum(B1)/N;
a = [a1; a];
b = [0 b];
v = [v vlr(1) 0];

end

