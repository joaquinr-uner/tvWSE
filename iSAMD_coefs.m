function vh = iSAMD_coefs(v,vnv,D,N,flag)
vh = zeros(1,2*sum(vnv));
vgam = v(end-D+1:end)./v(end-2*D+1:end-D);
for i=2:D
    if flag == 1
        A_l = v(i)*ones(1,vnv(i-1));
        gam_l = vgam(i);
    else
        A_l = v(end-2*D+i)*v(sum(vnv(1:i-2))+1:sum(vnv(1:i-1)));
        gam_l = vgam(i);
    end
    inter = floor(linspace(1,N,vnv(i-1)));
    t_l = inter(2:end-1);
    aux = [t_l,A_l,gam_l,i];
    vh(2*(sum(vnv(1:i-2)))+1:2*sum(vnv(1:i-1))) = aux;
end
end