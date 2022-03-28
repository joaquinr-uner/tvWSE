function [ti,alp,gam,e] = parse_coefs(N,v,D,vnv,t_corr)
ti = cell(1,D-1);
alp = cell(1,D-1);
gam = zeros(1,D-1);
e = zeros(1,D-1);
for i=1:D-1
    nv = vnv(i);
    vi = v(2*(sum(vnv(1:i-1)))+1:2*sum(vnv(1:i)));
    t_l = vi(1:vnv(i)-2);
    t_l = sort(round(t_l));
    if t_corr
        if t_l(1) <= 1
            t_l(1) = round((t_l(2)+1)/2);
        end
        
        if t_l(end) >= N
            t_l(end) = round((N+t_l(end-1))/2);
        end
        for j=2:length(t_l)-1
            if t_l(j)<=1 || t_l(j) >=N || t_l(j) == t_l(j+1)
                t_l(j) = round((t_l(j-1)+t_l(j+1))/2);
            end
        end
    end
    ti{i} = [1 t_l N];
    %ti{i} = t_l;
    alp{i} = vi(nv-1:2*nv-2);
    gam(i) = vi(2*nv-1);
    e(i) = vi(2*nv);
end

end