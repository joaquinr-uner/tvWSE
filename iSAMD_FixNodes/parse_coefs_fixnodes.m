function [t_ih,alp,gam,e] = parse_coefs_fixnodes(v,D,nv,N)

if length(nv) == 1
    vnv = nv*ones(1,D-1);
else
    vnv = nv;
end

alp = cell(1,D-1);
t_ih = cell(1,D-1);
gam = zeros(1,D-1);
e = zeros(1,D-1);
for i=1:D-1
    nv = vnv(i);
    vi = v(sum(vnv(1:i-1))+1+(i-1)*2:sum(vnv(1:i))+(i)*2);
    alp{i} = vi(1:nv);
    gam(i) = vi(nv+1);
    e(i) = vi(nv+2);
    t_ih{i} = linspace(1,N,nv);
end

end