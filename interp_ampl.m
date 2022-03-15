function A = interp_ampl(t,a,N,mInterp)
% Función que toma los coeficientes correspondientes del vector e interpola
% las amplitudes instantáneas utilizando splines cúbicos
% 
L = size(t,2);
A = zeros(L,N);

    for i=1:L
        x_l = t{i};
        y_l = a{i};
        A(i,:) = interp1(x_l,y_l,1:N,mInterp);
    end
end