function s = compute_splines(x,y,xq)
    [M,~] = size(x);
    N = length(xq);
    s = zeros(M,N);
    for i=1:M
       si = spline(x(i,:),y(i,:),xq);
       s(i,:) = si;
    end
end