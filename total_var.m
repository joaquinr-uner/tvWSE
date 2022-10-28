function TV = total_var(alpha_t)
    
    [M,N] = size(alpha_t);
    TV = zeros(1,M);

    for i=1:M
        TV(i) = sum(abs(diff(alpha_t(i,:))));
    end

end