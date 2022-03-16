function C = construct_dct(A,phi,D)
    [I,N] = size(A);
    
    C = zeros(N,2*sum(D));
    for i=1:I
       C(:,2*sum(D(1:i-1))+1:2*sum(D(1:i))) = cosenos(A(i,:),phi(i,:),D(i));
    end

end