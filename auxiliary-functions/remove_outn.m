function [ti2,alp2] = remove_outn(ti,alp,outn)
    D = length(ti);
    ti2 = cell(1,D);
    alp2 = cell(1,D);
    
    for i=1:D
        two = ti{i};
        alpwo = alp{i};
        tii = two(outn+1:end-outn);
        tii = tii - tii(1)+1;
        alpi = alpwo(outn+1:end-outn);
        
        ti2{i} = tii;
        alp2{i} = alpi;
    end 

end

