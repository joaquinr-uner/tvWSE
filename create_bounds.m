function [lb,ub] = create_bounds(vnv,vh,D,N)
    ub = zeros(1,2*sum(vnv));
    lb = zeros(1,2*sum(vnv));
    for i=1:D-1
       vhi = vh(2*(sum(vnv(1:i-1))+2*(i-1))+1:2*(sum(vnv(1:i))));
       ti = vhi(1:vnv(i)-2);
       %lbi = [1/N*ones(1,vnv(i)-2) zeros(1,vnv(i)) -Inf*ones(1,2)];
       lbi = [0 ti(1:end-1)+1/N -Inf*ones(1,vnv(i)+2)];
       %lbi = [ones(1,vnv(i)-2) zeros(1,vnv(i)) -Inf*ones(1,2)];
       ubi = [ti(2:end) 1 Inf*ones(1,vnv(i)+2)];
       
       lb(2*(sum(vnv(1:i-1)))+1:2*sum(vnv(1:i))) = lbi;
       ub(2*(sum(vnv(1:i-1)))+1:2*sum(vnv(1:i))) = ubi;
       
    end

end