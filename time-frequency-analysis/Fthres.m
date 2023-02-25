function F_thr = Fthres(F,mode)
% Perform thresholding of the STFT F. 
% Mode 1: Hard Thresholding
% Mode 2: Soft Thresholding
F_thr = zeros(size(F));

thr = 3*sqrt(2)*median(abs(real(F(:))))/0.6745;

[N,K] = size(F);

for i=1:N
     for j=1:K
         if (abs(F(i,j))>thr && mode==1)
             F_thr(i,j) = F(i,j);
         else if (abs(F(i,j))>thr && mode==2)
                 F_thr(i,j) = F(i,j) - thr;
             end
         end
     end
end


end
