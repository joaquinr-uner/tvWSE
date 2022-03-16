function A = interp_ampl_fixnodes(tn,v,N,vnv,D,mInterp)
% Función que toma los coeficientes correspondientes del vector e interpola
% las amplitudes instantáneas utilizando splines cúbicos
% 
A = zeros(D-1,N);

    for i=1:D-1
        nv = vnv(i);
        vi = v(sum(vnv(1:i-1))+1+(i-1)*2:sum(vnv(1:i))+(i)*2);
        y_l = vi(1:nv);
        A(i,:) = interp1(tn{i},y_l,1:N,mInterp);
        A(i,:) = A(i,:);
    end
end