function [alp,gam,e] = parse_coefs_fixnodes(v,D,vnv)
    alp = cell(1,D-1);
    gam = zeros(1,D-1);
    e = zeros(1,D-1);
    for i=1:D-1
        nv = vnv(i);
        vi = v(sum(vnv(1:i-1))+1+(i-1)*2:sum(vnv(1:i))+(i)*2);
        alp{i} = vi(1:nv);
        gam(i) = vi(nv+1);
        e(i) = vi(nv+2);
    end

end