function [lb,ub] = create_bounds(vnv,vh,D,N,outn)
%% Compute lower and upper bounds for the curve estimation algorithm.
% Inputs: 
%	  vnv: (D-1) integer vector with number of nodes for each harmonic component.
%	  vh: Initial coefficient vectors for tvWSE algorithm.
%	  D: number of harmonics of x.
% 	  N: length of signal.
% 	  outn: Number of outer nodes.
% Outputs: 
%	  lb: lower bound vector. 
%	  up: upper bound vector.

ub = zeros(1,2*sum(vnv)+4*outn*(D-1));
lb = zeros(1,2*sum(vnv)+4*outn*(D-1));
for i=1:D-1
    vhi = vh(2*(sum(vnv(1:i-1))+2*outn*(i-1))+1:2*(sum(vnv(1:i))+2*outn*(i)));
    ti = vhi(1:vnv(i)-2*(1-outn));
    
    lbi = [1/N, ti(1:end-1)+1/N -Inf*ones(1,vnv(i)+2*(outn+1))];
    ubi = [ti(2:end)-1/N, 1-1/N Inf*ones(1,vnv(i)+2*(outn+1))];

    if outn>0
        lbi(outn) = ti(outn);
        lbi(vnv(i)-2+2*outn) = ti(end-outn+1);
        ubi(outn) = ti(outn);
        ubi(vnv(i)-2+2*outn) = ti(end-outn+1);
    end

    %lbi = [zeros(1,outn-1) ti(1) ti(1:end-2)+1/N ti(end)*ones(1,outn) -Inf*ones(1,vnv(i)+2*outn+2)];
    %lbi = [ones(1,vnv(i)-2) zeros(1,vnv(i)) -Inf*ones(1,2)];
    
    lb(2*(sum(vnv(1:i-1))+2*outn*(i-1))+1:2*(sum(vnv(1:i))+2*outn*(i))) = lbi;
    ub(2*(sum(vnv(1:i-1))+2*outn*(i-1))+1:2*(sum(vnv(1:i))+2*outn*(i))) = ubi;
    
end

end