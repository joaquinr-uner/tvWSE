% This experiment implements that synthetic signal decomposition experiment
% detailed in Sec. 4.3 of "Fully Adaptive Time-Varying
% Wave-Shape Model: Applications in Biomedical Signal Processing".

addpath(genpath('/time-frequency-analysis'))
addpath(genpath('/auxiliary-functions'))

fs = 2000;
N = 2000;
fmax = 0.5;

t = 0:1/fs:N/fs-1/fs;
f = 0:fs/N:fs*fmax-fs/N;

A1 = sqrt(0.01*t)+1.1;
A2 = 2*log10((t+1.1))+0.5;

phi1 = 25*t + 5/(2*pi)*cos(2*pi*t);
phi2 = 100*t + 7*t.^2;

c1 = round(diff(phi1)*fs)+1;
c2 = round(diff(phi2)*fs)+1;
c1 = [c1(1), c1];
c2 = [c2(1), c2];g

D1 = 3;
D2 = 4;


sigma = 1e-4;
b = round(2.5/pi*sqrt(sigma/2)*N);

outn = 1;

alpha1 = zeros(D1,N);
alpha2 = zeros(D2,N); 

alpha1(1,:) = ones(1,N);
alpha1(2,:) = 0.5 + 0.25*cos(2*pi*3*t);
alpha1(3,:) = 0.3 + 0.2*exp(-(t-0.25).^2*10^2);

alpha2(1,:) = ones(1,N);
alpha2(2,:) = 0.5 + 0.4*t.^2;
alpha2(3,:) = 0.4 + 0.5*tanh(t-0.5);
alpha2(4,:) = 0.3 + 0.4*cos(2*pi*4*t);

e1 = [1 2.005, 3.003];
e2 = [1 2.002 3.002 3.998];

s1 = 0;

for i=1:D1
    s1 = s1 + alpha1(i,:).*cos(2*pi*e1(i)*phi1);
end
s1 = A1.*s1;

s2 = 0;

for i=1:D2
    s2 = s2 + alpha2(i,:).*cos(2*pi*e2(i)*phi2);
end

s2 = A2.*s2;

s_c = s1 + s2;

Nr = 200;
SNRs = [0,5,10,15,20];
for i=1:length(SNRs)

    SNR = SNRs(i);

    SNR1 = zeros(Nr,1);
    SNR1LR = SNR1;
    SNR1SAMD = SNR1;
    SNR1tvWSE = SNR1;
    
    SNR2 = SNR1;
    SNR2LR = SNR1;
    SNR2SAMD = SNR1;
    SNR2tvWSE = SNR1;

    Signals = zeros(Nr,N);
    
    S1LR = zeros(Nr,N);
    S2LR = zeros(Nr,N);
    
    S1SAMD = zeros(Nr,N);
    S2SAMD = zeros(Nr,N);
    
    S1tvWSE = zeros(Nr,N);
    S2tvWSE = zeros(Nr,N);
    
    vc = [6,8,10];
    parfor j=1:Nr
        r = 10^(-SNR/20)*std(s_c)*randn(size(s_c));

        snr1 = 20*log10(std(s1)/std(r));
        snr2 = 20*log10(std(s2)/std(r));
        
        s = s_c + r;

        [F, sF] = STFT_Gauss(s,N,1e-4,fmax);
        c1_est = ridge_ext(F,0.1,0.1,10,10);

        [A1_est,phi1_est] = extract_harmonics(F,sF,c1_est,b,b,1);

        r1_est = order_opt(s,40,A1_est,phi1_est,{'Wang'},vc);

        if r1_est<2
            r1_est = 2;
        end

        C1_est = construct_dct(A1_est,phi1_est,r1_est);
        v1_est = ((C1_est'*C1_est)\C1_est')*s';

        Np = 0.15*N;
        s_ext = extendSig(s,phi1_est,3,Np,'fw-bw');
        Next = length(s_ext);

        [Fext, sFext] = STFT_Gauss(s_ext,Next,sigma,fmax);

        c1ext = ridge_ext(Fext,0.1,0.1,10,10);
        [A1_ext,phi1_ext] = extract_harmonics(Fext,sFext,c1ext,b,b,1);

        C1_ext = construct_dct(A1_ext,phi1_ext,r1_est);
        v1_ext = ((C1_ext'*C1_ext)\C1_ext')*s_ext;

        s1n_ext = s_ext./(v1_ext(1)*A1_ext');
        C1n_ext = construct_dct(ones(1,Next),phi1_ext,r1_est);
        v1n_ext = ((C1n_ext'*C1n_ext)\C1n_ext')*s1n_ext;

        vnv1 = NNodes(0,r1_est,F,sF,c1_est,b,fs,0.8);

        vh1_ext = Init_tvWSE(v1n_ext,vnv1,r1_est,1,N,Next,outn);
        [lb1e,ub1e] = create_bounds(vnv1,vh1_ext,r1_est,Next,outn);

        [sn1_tvwse_ext,v_ie1_ext] = tvWSE(s1n_ext',ones(1,Next),phi1_ext,r1_est,vnv1,vh1_ext,'pchip',lb1e,ub1e,1,1,outn);
        s1_tvwse_ext = v1_ext(1)*A1_ext.*sn1_tvwse_ext;

        [ti,alpi,gamh,eh] = parse_coefs(Next,v_ie1_ext,r1_est,vnv1,0);


        sr_ext = s_ext-s1_tvwse_ext';
        [Fr_ext, sFr_ext] = STFT_Gauss(sr_ext,Next,sigma,fmax);

        c2ext = ridge_ext(Fr_ext,0.1,0.1,10,10);
        [A2_ext,phi2_ext] = extract_harmonics(Fr_ext,sFr_ext,c2ext,b,b,1);

        r2_est = order_opt(sr_ext,40,A2_ext,phi2_ext,{'Wang'},vc);
        
        if r2_est<2
            r2_est = 2;
        end
        
        C2_ext = construct_dct(A2_ext,phi2_ext,r2_est);
        v2_ext = ((C2_ext'*C2_ext)\C2_ext')*sr_ext;

        s2n_ext = sr_ext./(v2_ext(1)*A2_ext');
        C2n_ext = construct_dct(ones(1,Next),phi2_ext,r2_est);
        v2n_ext = ((C2n_ext'*C2n_ext)\C2n_ext')*s2n_ext;

        vnv2 = NNodes(0,r2_est,Fr_ext,sFr_ext,c2ext,b,fs,0.8);

        vh2_ext = Init_tvWSE(v2n_ext,vnv2,r2_est,1,N,Next,outn);
        [lb2e,ub2e] = create_bounds(vnv2,vh2_ext,r2_est,Next,outn);

        [sn2_tvwse_ext,v_ie2_ext] = tvWSE(s2n_ext',ones(1,Next),phi2_ext,r2_est,vnv2,vh2_ext,'pchip',lb2e,ub2e,1,1,outn);
        s2_tvwse_ext = v2_ext(1)*A2_ext.*sn2_tvwse_ext;

        s1_lr_ext = C1_ext*v1_ext;
        s2_lr_ext = C2_ext*v2_ext;

        [s_ext_samd,modes_ext_samd] = SAMD(s_ext',[A1_ext;A2_ext],[phi1_ext*2*pi;phi2_ext*2*pi],[r1_est,r2_est],1);

        [ti2,alpi2,gamh2,eh2] = parse_coefs(Next,v_ie2_ext,r2_est,vnv2,0);

        %%
        s1_lr = s1_lr_ext(Np+1:N+Np);
        s2_lr = s2_lr_ext(Np+1:N+Np);

        s1_samd = modes_ext_samd(1,Np+1:N+Np);
        s2_samd = modes_ext_samd(2,Np+1:N+Np);

        s1_tvwse = s1_tvwse_ext(Np+1:N+Np);
        s2_tvwse = s2_tvwse_ext(Np+1:N+Np);

        snr1ex_lr = 20*log10(std(s1)/std(s1-s1_lr'));
        snr2ex_lr = 20*log10(std(s2)/std(s2-s2_lr'));

        snr1ex_samd = 20*log10(std(s1)/std(s1-s1_samd));
        snr2ex_samd = 20*log10(std(s2)/std(s2-s2_samd));

        snr1ex_tvwse = 20*log10(std(s1)/std(s1-s1_tvwse));
        snr2ex_tvwse = 20*log10(std(s2)/std(s2-s2_tvwse));
        
        Signals(j,:) = s;

        SNR1(j) = snr1;
        SNR2(j) = snr2;

        S1LR(j,:) = s1_lr;
        S2LR(j,:) = s2_lr;        
        SNR1LR(j) = snr1ex_lr;
        SNR2LR(j) = snr2ex_lr;

        S1SAMD(j,:) = s1_samd;
        S2SAMD(j,:) = s2_samd;
        SNR1SAMD(j) = snr1ex_samd;
        SNR2SAMD(j) = snr2ex_samd;

        S1tvWSE(j,:) = s1_tvwse;
        S2tvWSE(j,:) = s2_tvwse;
        SNR1tvWSE(j) = snr1ex_tvwse;      
        SNR2tvWSE(j) = snr2ex_tvwse;
        
        fprintf('Denoising completed for signal 1. SNR = %.3f (LR), %.3f (SAMD), %.3f (tvWSE) \n',snr1ex_lr,snr1ex_samd,snr1ex_tvwse)
        fprintf('Denoising completed for signal 2. SNR = %.3f (LR), %.3f (SAMD), %.3f (tvWSE) \n',snr2ex_lr,snr2ex_samd,snr2ex_tvwse)
    end
    S = struct('Signals',Signals,'SNR1',SNR1,'SNR2',SNR2,'S1LR',S1LR,...
               'S2LR',S2LR,'SNR1LR',SNR1LR,'SNR2LR',SNR2LR,'S1SAMD',...
               S1SAMD,'S2SAMD',S2SAMD,'SNR1SAMD',SNR1SAMD,'SNR2SAMD',...
               SNR2SAMD,'S1tvWSE',S1tvWSE,'S2tvWSE',S2tvWSE,'SNR1tvWSE',...
               SNR1tvWSE,'SNR2tvWSE',SNR2tvWSE);
    save(['Results_NoiseMultiComp_' num2str(SNR) 'dB'],'S')
end