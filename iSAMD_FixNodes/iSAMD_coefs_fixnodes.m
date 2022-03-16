function [vh,tn] = iSAMD_coefs_fixnodes(v,vnv,D,N,flag)

vh = zeros(1,sum(vnv)+2*(D-1));
tn = cell(1,D-1);
for i=2:D
    if flag == 1
        A_l = v(i)*ones(1,vnv(i-1));
        gam_l = v(i)/v(D+i);
    else
        ah = v(sum(vnv(1:i-2))+1:sum(vnv(1:i-1)));
        vlr = v(end-2*D+1:end);
        A_l = vlr(i)*v(sum(vnv(1:i-2))+1:sum(vnv(1:i-1)));
        gam_l = vlr(i+D)/vlr(i);
    end
    inter = floor(linspace(1,N,vnv(i-1)));
    tn{i-1} = inter;
    aux = [A_l,gam_l,i];
    vh(sum(vnv(1:i-2))+1+(i-2)*2:sum(vnv(1:i-1))+(i-1)*2) = aux;
end