function vh = Init_tvWSE(v,vnv,D,flag,N,Ne,t_init,outn)
% Constructs the initial coefficient vector to be used in the adaptive 
% time-varying wave-shape extraction (tvWSE) algorithm. See [1] for details.
% Inputs:
%	  - v: WSF linear regression coefficients for signal x. 
% 	  - vnv: (D-1) integer vector with number of nodes for each harmonic component.
%	  - D: number of harmonics of x.
%	  - flag: 
%	  - N: length of signal.
%	  - Ne: length of the extended signal.
%	  - t_init: cell structure with initial node locations. If equal to 0, 
%		    node locations are computed normally.
%	  - outn: Number of outer nodes.
% Outputs:
%	  - vh: Initial coefficient vectors for tvWSE algorithm.

if nargin<4
    flag = 1;
end
if nargin<5
    N = 1;
    Ne = 1;
end
if nargin<7
    t_init = 0;
end
if nargin<8 && Ne/N>1
    outn = 1;
else
    if nargin<8 && Ne/N==1
        outn = 0;
    end
end

r_ext = Ne/N;
vh = zeros(1,2*sum(vnv)+4*outn*(D-1));
vgam = v(end-D+1:end)./v(end-2*D+1:end-D);
for i=2:D
    if flag == 1
        A_l = v(i)*ones(1,vnv(i-1)+2*outn);
        gam_l = vgam(i);
    else
        A_l = v(sum(vnv(1:i-2)+2*outn)+1:sum(vnv(1:i-1)+2*outn));
        %A_l = v(sum(vnv(1:i-2))+1:sum(vnv(1:i-1)));
        gam_l = vgam(i);
        %gam_l = 0;
    end
    if iscell(t_init)
        inter = t_init{i-1};
    else
        %inter = floor(linspace(1,N,vnv(i-1)));
        if outn>0
            ti = 0.5*(1-1/r_ext+1/Ne);
            te = 1-0.5*(1-1/r_ext);
            outni = linspace(0,ti,outn+1);
            outne = linspace(te,1,outn+1);
            t_in = linspace(ti,te,vnv(i-1));
            inter = [outni(2:end-1) t_in outne(2:end-1)];
        %inter = linspace(0,1,vnv(i-1));
        else
            inter = linspace(0,1,vnv(i-1));
            inter = inter(2:end-1);
    end
    t_l = inter;
    aux = [t_l,A_l,gam_l,i];
    vh(2*(sum(vnv(1:i-2))+2*outn*(i-2))+1:2*(sum(vnv(1:i-1))+2*outn*(i-1))) = aux;
end
end