function [ti,alp,gam,e] = parse_coefs(N,v,D,nv,t_corr,outn)
%% Parse Optimal Coefficient vector of the adaptive
% time-varying wave-shape extraction (tvWSE) algorithm.
% Inputs:
%	  - N: length of signal.
%	  - v: Optimal Coefficient Vector from the tvWSE algorithm.
%	  - D: number of harmonics of x.
%	  - nv: (D-1) integer vector with number of nodes for each harmonic component.
%	  - t_corr: 
%	  - outn: Number of outer nodes.
% Outputs:
%	  - ti: cell structure with node locations for optimal nodes.
%	  - alp: cell structure with node amplitudes for optimal nodes.
%	  - gam: optimal phase-shift coefficients gamma_l.
%	  - e: optimal non-integer phase coefficients e_l.
if length(nv) == 1
    vnv = nv*ones(1,D-1);
else
    vnv = nv;
end

ti = cell(1,D-1);
alp = cell(1,D-1);
gam = zeros(1,D-1);
e = zeros(1,D-1);
for i=1:D-1
    nv = vnv(i);
    vi = v(2*sum(vnv(1:i-1))+4*outn*(i-1)+1:2*sum(vnv(1:i))+4*outn*i);
    t_l = vi(1:vnv(i)-2+2*outn);
    [t_l, indx_l] = sort(t_l);
    %t_l = round(t_l);
    if t_corr
        if t_l(1) <= 1
            t_l(1) = round((t_l(2)+1)/2);
        end
        
        if t_l(end) >= N
            t_l(end) = round((N+t_l(end-1))/2);
        end
        for j=2:length(t_l)-1
            if t_l(j)<=1 || t_l(j) >=N || t_l(j) == t_l(j+1) || t_l(j) == t_l(j-1)
                t_l(j) = round((t_l(j-1)+t_l(j+1))/2);
            end
        end
    end
    if t_l(end)<N/2
        ti{i} = [1 N*t_l N];
    else
        ti{i} = [1 t_l N];
    end
    %ti{i} = t_l;
    aux = vi(nv-2+2*outn+1:2*nv-2+4*outn);
    aux_in = aux(2:end-1);
    aux(2:end-1) = aux_in(indx_l);
    alp{i} = aux;
    gam(i) = vi(end-1);
    e(i) = vi(end);
end

end