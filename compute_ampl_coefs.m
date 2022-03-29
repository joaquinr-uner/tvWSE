function [a, b, alp, v] = compute_ampl_coefs(ti_h,Al,gamh,eh,N,D,nv)

if length(nv) == 1
    vnv = nv*ones(1,D-1);
else
    vnv = nv;
end

alp = zeros(D-1,N);
a = sum(Al,2)/N;
b = gamh.*a';
v = [];
for i=1:D-1
    alp(i,:) = Al(i,:)/a(i);
    aux = [ti_h{i} alp(i,ti_h{i}) a(i) b(i) eh(i)];
    v = [v aux];
end


end

