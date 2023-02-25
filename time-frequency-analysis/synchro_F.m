function [SF] = synchro_F(F,fmax)
% synchro_F applies the synchrosqueezing operator on the STFT


[K, N] = size(F);
op_frec = -1/(2*pi)*diff(unwrap(angle(F')));
op_frec = [op_frec; op_frec(end,:)]; op_frec = op_frec';


SF = zeros(size(F));

for k=1:K
    for n=1:N
        aux = round(K/fmax*op_frec(k,n));
        if(aux>=1 && aux<=K)
            SF(aux,n) = SF(aux,n) + F(k,n);
        end
    end
end