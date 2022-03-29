function [a, b, alp, ti_h, eh] = compute_ampl_coefs(v,N,D,vnv,mInterp)

[ti_h,alph,gamh,eh] = parse_coefs(N,v,D,vnv,1);
Al = interp_ampl(ti_h,alph,N,mInterp);

alp = zeros(D-1,N);
a = sum(Al,2)/N;
b = gamh.*a;
for i=1:D-1
    alp(i,:) = Al(i,:)/a(i);
end


end

