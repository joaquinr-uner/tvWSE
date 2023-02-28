% This experiment implements that synthetic signal denoising experiment
% detailed in Sec. 4.2 of "Fully Adaptive Time-Varying
% Wave-Shape Model: Applications in Biomedical Signal Processing".

addpath(genpath('time-frequency-analysis'))
addpath(genpath('auxiliary-functions'))
mInterp = 'pchip';

fs = 2000;
N = 2000;
T = N/fs;
t = 0:1/fs:T-1/fs;
f = 0:fs/N:fs/2-fs/N;

D = 3;
alpha = zeros(D,N);

A = 0.1*sqrt(t+1);
phi = 40*t+5/(2*pi)*sin(2*pi*t);

alpha(1,:) = ones(1,N);

B2 = zeros(4,N);
B3 = zeros(4,N);

B2(1,:) = 0.5 + 0.25*cos(2*pi*3*t);
B2(2,:) = 0.5 + 0.3*tanh(25*(t-0.5));
B2(3,:) = 0.5*linspace(0,T,N)+0.1;
B2(4,:) = 0.5*linspace(0,T,N)+0.1;

B3(1,:) = 0.3 + 0.25*cos(2*pi*4*t);
B3(2,:) = 0.4*t + 0.4*exp(-(100*(t-0.25).^2));
B3(3,:) = 0.3 + 0.2*cos(2*pi*4*t);
B3(4,:) = 0.3 + 0.2*tanh(25*(t-0.75));

B2labs = {'Cosine','Tanh','Linear','Linear'};
B3labs = {'Cosine','Bump','Cosine','Tanh'};


e = [1, 2.005, 2.995];

Crit = {'Wang'};
Cparams = struct('c',[6,8,10,12]);
s_c = 0;
Nr = 100;
SNRs = [20,15,10,5,0];
Np = 0.1*N;
sigma = 1e-4;
fmax = 0.5;
b = round(3/pi*sqrt(sigma/2)*N);

for i=1:size(B2,1)
        
        alpha(2,:) = B2(i,:);
        alpha(3,:) = B3(i,:);
        s_c = alpha(1,:).*cos(2*pi*phi);
        for l=2:D
            s_c = s_c + alpha(l,:).*(cos(2*pi*e(l)*phi));
        end
        s_c = A.*s_c;
        
        
        E_lr = zeros(Nr,1);
        E_samd = zeros(Nr,1);
        E_tvWSE = zeros(Nr,1);
        E_ht = zeros(Nr,1);
        E_st = zeros(Nr,1);
        SNR_lr = zeros(Nr,1);
        SNR_samd = zeros(Nr,1);
        SNR_tvWSE = zeros(Nr,1);
        SNR_ht = zeros(Nr,1);
        SNR_st = zeros(Nr,1);
        T_samd = zeros(Nr,1);
        T_tvWSE = zeros(Nr,1);
        Signals = zeros(Nr,N);
        S_tvWSE = zeros(Nr,N);
        S_LR = zeros(Nr,N);
        S_SAMD = zeros(Nr,N);
        S_Hard = zeros(Nr,N);
        S_Soft = zeros(Nr,N);
        alp2est = zeros(Nr,N);
        alp3est = zeros(Nr,N);
        B2est = zeros(Nr,N);
        B3est = zeros(Nr,N);
        Ropt = zeros(Nr,1);
                
        for m=1:length(SNRs)
            SNR = SNRs(m);
            fprintf(['Denoising of time-varying waveshapes with B2 = ' B2labs{i} ' and B3 = ' B3labs{i} '. SNR_in = ' num2str(SNR) '...\n'])
            for k=1:Nr
                
                r = 10^(-SNR/20)*std(s_c)*randn(size(s_c));
                
                s = s_c + r;
                
                [F, sF] = STFT_Gauss(s,N,sigma,fmax);
                c = ridge_ext(F,0.1,0.1,10,10);
                
                r_max = floor(0.5*N/max(c));
                
                
                [A_est, phi_est] = extract_harmonics(F,sF,c,b,b,1);
                
                r_opt = order_opt(s',r_max,A_est,phi_est,Crit,Cparams,F);
                                
                if r_opt<2
                    r_opt = 2;
                end
                C = makeC(A_est,phi_est,r_opt);
                v = ((C'*C)\C')*s';
                sn = s./(v(1)*A_est);
            
                [Fn,sFn] = STFT_Gauss(sn,N,sigma,fmax);
            
                nv = 0;


                vnv = NNodes(nv,r_opt,Fn,sFn,c,b,fs,0.9);      
                
                                                
                cycls = 3;
                s_ext = extendSig(s,phi_est,cycls,Np,'fw-bw');        
                Next = N+2*Np;
                text = 0:1/fs:Next/fs-1/fs;
                fext = 0:fs/Next:fs*fmax-fs/Next;
                
                
                [Fext,sFext] = STFT_Gauss(s_ext,Next,sigma,fmax);
                cext = ridge_ext(Fext,0.1,0.1,10,10);
                be = round(3/pi*sqrt(sigma/2)*Next);
                [A_ext,phi_ext] = extract_harmonics(Fext,sFext,cext,be,be,1);
                                
                C = construct_dct(A_ext,phi_ext,r_opt);
                ve = ((C'*C)\C')*s_ext;
                se_lr = C*ve;
                se_lr = se_lr';
                tic;
                
                [se_samd,~,v_e] = my_SAMD(s_ext,A_ext,phi_ext*(2*pi),r_opt,1);
                t_samd = toc;
                                
                B1e = ve(1)*A_ext;
                se_norm = s_ext'./B1e;
                se_norm = se_norm - mean(se_norm);
                
                Cn = construct_dct(ones(1,Next),phi_ext,r_opt);
                vn = ((Cn'*Cn)\Cn')*se_norm';
                
                
                [Fen,sFen] = STFT_Gauss(se_norm,Next,sigma,fmax);
                
                vh = Init_tvWSE(vn,vnv,r_opt,1,N,Next);
                [lb,ub] = create_bounds(vnv,vh,r_opt,Next,1);
                options = optimoptions('lsqcurvefit','Algorithm','levenberg-marquardt',...
                'MaxFunctionEvaluations',400*length(vh),'Display','off');

                tic;
                [se_tvwsen, v_ie] = tvWSE(se_norm,ones(1,Next),phi_ext,r_opt,vnv,vh,mInterp,lb,ub,1,1,1,options);
                
                se_tvwse = B1e.*se_tvwsen;
                
                t_tvwse = toc;
                
                F_hthr = Fthres(Fext,1);
                F_sthr = Fthres(Fext,2);
                se_hard= 2/max(sFext)*real(sum(F_hthr,1));
                se_soft = 2/max(sFext)*real(sum(F_sthr,1));
                se_soft = se_soft - mean(se_soft);
                
                s_lr = se_lr(Np+1:Next-Np);
                s_samd = se_samd(Np+1:Next-Np);
                s_tvwse = se_tvwse(Np+1:Next-Np);
                s_hard = se_hard(Np+1:Next-Np);
                s_soft = se_soft(Np+1:Next-Np);
                                                
                e_lr = norm(s_c - s_lr)/norm(s_c);
                e_samd = norm(s_c - s_samd)/norm(s_c);
                e_tvwse = norm(s_c - s_tvwse)/norm(s_c);
                e_ht = norm(s_c - s_hard)/norm(s_c);
                e_st = norm(s_c - s_soft)/norm(s_c);
                
                snr_lr = 20*log10(std(s_c)/std(s_c-s_lr));
                snr_samd = 20*log10(std(s_c)/std(s_c-s_samd));
                snr_tvwse = 20*log10(std(s_c)/std(s_c-s_tvwse));
                snr_ht = 20*log10(std(s_c)/std(s_c-s_hard));
                snr_st = 20*log10(std(s_c)/std(s_c-s_soft));
                
                E_lr(k,:) = e_lr;
                E_samd(k,:) = e_samd;
                E_tvWSE(k,:) = e_tvwse;
                E_ht(k,:) = e_ht;
                E_st(k,:) = e_st;
                SNR_lr(k,:) = snr_lr;
                SNR_samd(k,:) = snr_samd;
                SNR_tvWSE(k,:) = snr_tvwse;
                SNR_ht(k,:) = snr_ht;
                SNR_st(k,:) = snr_st;
                                
                T_samd(k) = t_samd;
                T_tvWSE(k) = t_tvwse;
                Signals(k,:) = s;
                S_tvWSE(k,:) = s_tvwse;
                S_LR(k,:) = s_lr;
                S_SAMD(k,:) = s_samd;
                S_Hard(k,:) = s_hard;
                S_Soft(k,:) = s_soft;
                Ropt(k,:) = r_opt;
                
                fprintf(['SNR_out (LR) = ' num2str(snr_lr) '. SNR_out (SAMD) = ' num2str(snr_samd) ...
                    '. SNR_out (tvWSE) = ' num2str(snr_tvwse) ...
                    '. SNR_out (HT) = ' num2str(snr_ht) ', SNR_out (ST) = ' num2str(snr_st) '...\n'])
            end
            
            S = struct('CleanS',s_c,'Signals',Signals,'S_LR',S_LR,'S_SAMD',S_SAMD,'S_tvWSE',...
                S_tvWSE,'S_Hard',S_Hard,'S_Soft',S_Soft,'Errors_tvWSE',...
                E_tvWSE,'Errors_SAMD',E_samd,'Errors_LR',E_lr,'Errors_HT',...
                E_ht,'Errors_ST',E_st,'SNR_LR',SNR_lr,'SNR_SAMD',SNR_samd,...
                'SNR_tvWSE',SNR_tvWSE,'SNR_HT',SNR_ht,'SNR_ST',SNR_st,...
                'Times_SAMD',T_samd,'Times_tvWSE',T_tvWSE);
            save(fullfile(drt_r,['Results_DenoisingExt_B2_' B2labs{i} '_B3_' B3labs{i} '_' num2str(SNR) 'dB_16Nov22.mat']),'S')
        end
end
