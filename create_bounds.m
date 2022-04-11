function [lb,ub] = create_bounds(vnv,D,N)
    ub = zeros(1,2*sum(vnv));
    lb = zeros(1,2*sum(vnv));
    for i=1:D-1
       lbi = [ones(1,vnv(i)-2) -Inf*ones(1,vnv(i)+2)];
       %lbi = [ones(1,vnv(i)-2) zeros(1,vnv(i)) -Inf*ones(1,2)];
       ubi = [N*ones(1,vnv(i)-2) Inf*ones(1,vnv(i)+2)];
       
       lb(2*(sum(vnv(1:i-1)))+1:2*sum(vnv(1:i))) = lbi;
       ub(2*(sum(vnv(1:i-1)))+1:2*sum(vnv(1:i))) = ubi;
       
    end

end