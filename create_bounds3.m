function [lb,ub] = create_bounds3(vnv,vh,D,N,outn)
ub = zeros(1,2*sum(vnv)+4*outn*(D-1));
lb = zeros(1,2*sum(vnv)+4*outn*(D-1));
for i=1:D-1
    vhi = vh(2*(sum(vnv(1:i-1))+2*outn*(i-1))+1:2*(sum(vnv(1:i))+2*outn*(i)));
    ti = vhi(outn:vnv(i)+outn-1);
    
    %lbi = [1/N*ones(1,vnv(i)-2) zeros(1,vnv(i)) -Inf*ones(1,2)];
    lbi = [zeros(1,outn-1) ti(1) ti(1:end-2)+1/N ti(end)*ones(1,outn) -Inf*ones(1,vnv(i)+2*outn+2)];
    %lbi = [ones(1,vnv(i)-2) zeros(1,vnv(i)) -Inf*ones(1,2)];
    ubi = [ti(1)*ones(1,outn) ti(3:end)-1/N ti(end) ones(1,outn-1) Inf*ones(1,vnv(i)+2*outn+2)];
    
    lb(2*(sum(vnv(1:i-1))+2*outn*(i-1))+1:2*(sum(vnv(1:i))+2*outn*(i))) = lbi;
    ub(2*(sum(vnv(1:i-1))+2*outn*(i-1))+1:2*(sum(vnv(1:i))+2*outn*(i))) = ubi;
    
end

end