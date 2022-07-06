function vh = iSAMD_Coefs2(v,vnv,D,flag,N,Next,t_init)
% Constructs the coefficient vector to be used in the iSAMD algorithm
% with the extended signal
if nargin<4
    flag = 1;
else
    if nargin<5
        N = 1;
        Next = 1;
    else
        if nargin<7
            t_init = 0;
        end
    end
end
vh = zeros(1,2*sum(vnv)+4*(D-1));
vgam = v(end-D+1:end)./v(end-2*D+1:end-D);
for i=2:D
    if flag == 1
        A_l = v(i)*ones(1,vnv(i-1)+2);
        gam_l = vgam(i);
    else
        A_l = v(sum(vnv(1:i-2)+2)+1:sum(vnv(1:i-1)+2));
        %A_l = v(sum(vnv(1:i-2))+1:sum(vnv(1:i-1)));
        gam_l = vgam(i);
        %gam_l = 0;
    end
    if iscell(t_init)
        inter = t_init{i-1};
    else
        %inter = floor(linspace(1,N,vnv(i-1)));
        ti = ((Next-N+1)/2)/Next;
        te = 1-((Next-N)/2)/Next;
        inter = linspace(ti,te,vnv(i-1));
        %inter = linspace(0,1,vnv(i-1));
    end
    t_l = inter;
    aux = [t_l,A_l,gam_l,i];
    vh(2*(sum(vnv(1:i-2)+(i-1)))+1:2*sum(vnv(1:i-1))+4*(i-1)) = aux;
end
end
