function A = extract_ridge(F,N,c)
    [~,N] = size(F);
    A = zeros(1,N);
    for i=1:N
        A(i) = 1/N*abs(F(c(i),i)).^2;
    end
end