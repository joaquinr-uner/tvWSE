% This experiment implements that synthetic signal segmentation experiment
% detailed in Sec. 4.4 of "Fully Adaptive Time-Varying
% Wave-Shape Model: Applications in Biomedical Signal Processing".
addpath(genpath('/time-frequency-analysis'))
addpath(genpath('/auxiliary-functions'))

mInterp = 'pchip';


fs = 2000;
N = 2000;
T = N/fs;
t = 0:1/fs:T-1/fs;
f = 0:fs/N:fs/2-fs/N;

A = 0.1*sqrt(t+1);
phi = 40*t+5/(2*pi)*sin(2*pi*t);


e = [1, 2.005, 2.995, 4.005, 4.995, 5.995];

Crit = {'Wang'};
Cparams = struct('c',[6,8,10,12]);
Nr = 100;
SNRs = [20,15,10,5,0];
Np = 0.1*N;
sigma = 1e-4;
fmax = 0.5;
b = round(3/pi*sqrt(sigma/2)*N);

for m=1:length(SNRs)
    SNR = SNRs(m);

    E_lr = zeros(Nr,2);
    E_samd = zeros(Nr,2);
    E_tvWSE = zeros(Nr,2);
    E_ht = zeros(Nr,2);
    E_st = zeros(Nr,2);
    SNR_lr = zeros(Nr,2);
    SNR_samd = zeros(Nr,2);
    SNR_tvWSE = zeros(Nr,2);
    SNR_ht = zeros(Nr,2);
    SNR_st = zeros(Nr,2);

    T_samd = zeros(Nr,1);
    T_tvWSE = zeros(Nr,1);
    VD = zeros(Nr,1);
    Ropt = zeros(Nr,1);
    Tt = zeros(Nr,1);
    S_C = zeros(Nr,N);
    Signals = zeros(Nr,N);
    S_tvWSE = zeros(Nr,N);
    S_LR = zeros(Nr,N);
    S_SAMD = zeros(Nr,N);
    S_Hard = zeros(Nr,N);
    S_Soft = zeros(Nr,N);
    Alp = zeros(Nr,6,N);
    AlpEst = zeros(Nr,6,N);
    ACoefs = zeros(Nr,6);
    Nodes = cell(Nr,1);
    fprintf(['Segmentation of time-varying waveshapes with SNR_in = ' num2str(SNR) '...\n'])
    parfor k=1:Nr
        D = 2 + randi(4,[1,1]);
        t_t = round(0.1*N + randi(0.8*N,[1,1]));
        alpha = zeros(6,N);
        alpha(1,:) = ones(1,N);
        
        for i=2:D
            alpha(i,:) = (0.1 + 0.4*rand(1,1)) + (0.1 + 0.25*rand(1,1))*tanh(50*(t-t(t_t)));
        end

        s_c = alpha(1,:).*cos(2*pi*phi);
        for l=2:D
            s_c = s_c + alpha(l,:).*(cos(2*pi*e(l)*phi));
        end
        s_c = A.*s_c;

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

        [ti,ql,gamh,eh] = parse_coefs(Next,v_ie,r_opt,vnv,0,1);
        [a,~,Al,alp] = compute_hafs(ti,ql,gamh,'pchip',1,B1e);
        
        alpha_est = ones(6,N);
        VA = zeros(1,6);
        VA(1:r_opt) = a;

        for i=2:r_opt
            alpha_est(i,:) = Al(i-1,Np+1:N+Np);
        end

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

        % Errors Pre-Transition
        epr_lr = norm(s_c(1:t_t) - s_lr(1:t_t))/norm(s_c(1:t_t));
        epr_samd = norm(s_c(1:t_t) - s_samd(1:t_t))/norm(s_c(1:t_t));
        epr_tvwse = norm(s_c(1:t_t) - s_tvwse(1:t_t))/norm(s_c(1:t_t));
        epr_ht = norm(s_c(1:t_t) - s_hard(1:t_t))/norm(s_c(1:t_t));
        epr_st = norm(s_c(1:t_t) - s_soft(1:t_t))/norm(s_c(1:t_t));

        % SNR Pre-Transition
        snrpr_lr = 20*log10(std(s_c(1:t_t))/std(s_c(1:t_t)-s_lr(1:t_t)));
        snrpr_samd = 20*log10(std(s_c(1:t_t))/std(s_c(1:t_t)-s_samd(1:t_t)));
        snrpr_tvwse = 20*log10(std(s_c(1:t_t))/std(s_c(1:t_t)-s_tvwse(1:t_t)));
        snrpr_ht = 20*log10(std(s_c(1:t_t))/std(s_c(1:t_t)-s_hard(1:t_t)));
        snrpr_st = 20*log10(std(s_c(1:t_t))/std(s_c(1:t_t)-s_soft(1:t_t)));

        %Error Post-Transition
        eps_lr = norm(s_c(t_t+1:end) - s_lr(t_t+1:end))/norm(s_c(t_t+1:end));
        eps_samd = norm(s_c(t_t+1:end) - s_samd(t_t+1:end))/norm(s_c(t_t+1:end));
        eps_tvwse = norm(s_c(t_t+1:end) - s_tvwse(t_t+1:end))/norm(s_c(t_t+1:end));
        eps_ht = norm(s_c(t_t+1:end) - s_hard(t_t+1:end))/norm(s_c(t_t+1:end));
        eps_st = norm(s_c(t_t+1:end) - s_soft(t_t+1:end))/norm(s_c(t_t+1:end));

        %SNR Post-Transition
        snrps_lr = 20*log10(std(s_c(t_t+1:end))/std(s_c(t_t+1:end)-s_lr(t_t+1:end)));
        snrps_samd = 20*log10(std(s_c(t_t+1:end))/std(s_c(t_t+1:end)-s_samd(t_t+1:end)));
        snrps_tvwse = 20*log10(std(s_c(t_t+1:end))/std(s_c(t_t+1:end)-s_tvwse(t_t+1:end)));
        snrps_ht = 20*log10(std(s_c(t_t+1:end))/std(s_c(t_t+1:end)-s_hard(t_t+1:end)));
        snrps_st = 20*log10(std(s_c(t_t+1:end))/std(s_c(t_t+1:end)-s_soft(t_t+1:end)));

        E_lr(k,:) = [epr_lr,eps_lr];
        E_samd(k,:) = [epr_samd,eps_samd];
        E_tvWSE(k,:) = [epr_tvwse,eps_tvwse];
        E_ht(k,:) = [epr_ht,eps_ht];
        E_st(k,:) = [epr_st,eps_st];
        SNR_lr(k,:) = [snrpr_lr,snrps_lr];
        SNR_samd(k,:) = [snrpr_samd,snrps_samd];
        SNR_tvwse(k,:) = [snrpr_tvwse,snrps_tvwse];
        SNR_ht(k,:) = [snrpr_ht,snrps_ht];
        SNR_st(k,:) = [snrpr_st,snrps_st];

        T_samd(k) = t_samd;
        T_tvWSE(k) = t_tvwse;
        Signals(k,:) = s;
        S_C(k,:) = s_c;
        S_tvWSE(k,:) = s_tvwse;
        S_LR(k,:) = s_lr;
        S_SAMD(k,:) = s_samd;
        S_Hard(k,:) = s_hard;
        S_Soft(k,:) = s_soft;
        Alp(k,:,:) = alpha;
        AlpEst(k,:,:) = alpha_est;
        ACoefs(k,:) = VA;
        Ropt(k) = r_opt;
        VD(k) = D;
        Tt(k) = t_t;
        Nodes{k} = v_ie;
        fprintf(['SNR_out (LR) = ' num2str(snrpr_lr) '. SNR_out (SAMD) = ' num2str(snrpr_samd) ...
            '. SNR_out (tvWSE) = ' num2str(snrpr_tvwse) ...
            '. SNR_out (HT) = ' num2str(snrpr_ht) ', SNR_out (ST) = ' num2str(snrpr_st) '...\n'])
        fprintf(['SNR_out (LR) = ' num2str(snrps_lr) '. SNR_out (SAMD) = ' num2str(snrps_samd) ...
            '. SNR_out (tvWSE) = ' num2str(snrps_tvwse) ...
            '. SNR_out (HT) = ' num2str(snrps_ht) ', SNR_out (ST) = ' num2str(snrps_st) '...\n'])
    end

    S = struct('CleanS',S_C,'Signals',Signals,'S_LR',S_LR,'S_SAMD',S_SAMD,'S_tvWSE',...
        S_tvWSE,'S_Hard',S_Hard,'S_Soft',S_Soft,'Errors_tvWSE',...
        E_tvWSE,'Errors_SAMD',E_samd,'Errors_LR',E_lr,'Errors_HT',...
        E_ht,'Errors_ST',E_st,'SNR_LR',SNR_lr,'SNR_SAMD',SNR_samd,...
        'SNR_tvWSE',SNR_tvWSE,'SNR_HT',SNR_ht,'SNR_ST',SNR_st,...
        'Times_SAMD',T_samd,'Times_tvWSE',T_tvWSE,'ALP',Alp,...
        'AlpEst',AlpEst,'ACoefs',ACoefs,'R_opt',Ropt,'D',VD,'Tt',Tt);
    save(fullfile(drt_r,['Results_Segmentation_' num2str(SNR) 'dB.mat']),'S')
    save(fullfile(drt_r,['Nodes_Segmentation_' num2str(SNR) 'dB.mat']),'Nodes')
end