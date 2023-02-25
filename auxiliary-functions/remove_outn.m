function [ti2,alp2] = remove_outn(ti,alp,outn)
% Remove outer nodes from node location (ti) and amplitude (alp) vectors.
% Inputs:
% 	 - ti: cell structure with the node locations t_{i,l} for each harmonic.
%	 - alp: cell structure with the node amplitudes alp_{i,l} for each harmonic.
% 	 - outn: number of outer nodes at the start and end of the signal.
% Outputs:
%	 - ti2: cell structure with the node locations t_{i,l} for each harmonic, with outer nodes removed.
%	 - alp2: cell structure with the node amplitudes alp_{i,l} for each harmonic, with outer nodes removed.

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

