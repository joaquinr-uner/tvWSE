function [lb,ub] = create_bounds2(vnv,vh,D,N)
ub = zeros(1,2*sum(vnv)+4*(D-1));
lb = zeros(1,2*sum(vnv)+4*(D-1));
for i=1:D-1
    vhi = vh(2*(sum(vnv(1:i-1))+2*(i-1))+1:2*(sum(vnv(1:i))+2*i));
    t1i = vhi(1);
    t1f = vhi(vnv(i));
    %lbi = [1/N*ones(1,vnv(i)-2) zeros(1,vnv(i)) -Inf*ones(1,2)];
    lbi = [t1i t1i+1/N*ones(1,vnv(i)-2) t1f -Inf*ones(1,vnv(i)+4)];
    %lbi = [ones(1,vnv(i)-2) zeros(1,vnv(i)) -Inf*ones(1,2)];
    ubi = [t1i (t1f-1/N)*ones(1,vnv(i)-2) t1f Inf*ones(1,vnv(i)+4)];
    
    lb(2*(sum(vnv(1:i-1))+2*(i-1))+1:2*(sum(vnv(1:i))+2*i)) = lbi;
    ub(2*(sum(vnv(1:i-1))+2*(i-1))+1:2*(sum(vnv(1:i))+2*i)) = ubi;
    
end

end