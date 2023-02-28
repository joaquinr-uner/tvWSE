% This Script generates Fig. 7 from the paper "Fully Adaptive Time-Varying
% Wave-Shape Model: Applications in Biomedical Signal Processing"

addpath(genpath('auxiliary-functions'))
addpath(genpath('time-frequency-analysis'))
[x, hdr, label, fs, scle, offs] = read_edf('eeg44.edf');

fs = 256;
index = 343*fs+1:374*fs;
x = x{16}(index) - x{10}(index); 
x = double(x);
x = x - mean(x);
N = length(x);
t = 0:1/fs:N/fs-1/fs;
%%

% Parameters for STFT
bt = 1:N;
redun = 8;
gamma = 0;
sigma = 0.1;
fmax = 0.05;
ft = 1/redun:1/redun:round(fmax*N);

[STFT,SST1,SST2] = sstn_test_modL_new(x,gamma,sigma,ft,bt,redun);

RTF = SST2;
d = redun;
jump = redun/2; 
c1 = exridge(RTF,0,0,jump);
aux = zeros(size(t));
for j = 1:N
    aux(j) = sum(RTF(c1(j)-d:c1(j)+d,j));
    RTF(c1(j)-d:c1(j)+d,j) = 0;
end;
A = abs(aux); 
phi = phase(aux)/(2*pi);

D1 = 6; 

tic;
C = [makeC(A,phi,D1)];

C1_aux = makeC(ones(1,1000),phi(5*fs+1):1/1000:phi(5*fs+1)+1-1/1000,D1);

coef_est = ((C'*C)\C')*x';
f_est_LR = C*coef_est;
f1_est_LR = C(:,1:2*D1)*coef_est(1:2*D1);
tiempo_LR = toc

WSF1 = C1_aux*coef_est(1:2*D1);

order_Phi = 1; 
warning('off')
fprintf('Computing SAMD')
tic;
[s_samd,b_e] = one_wsf_nonlinear(x,A,phi*2*pi,D1,order_Phi);
tiempo_NLR = toc
warning('on')

figure(1)
subplot(311)
hold on
plot(t,s_samd)
hold off
%%
opt.maxiter = 10000;
opt.eps_error = 1e-12;
opt.show = 0;
opt.nknots = 20;
opt.knotremoval_factor= 1.01;
opt.order = 3;
opt.eps_diff = opt.eps_error;
    
ins_pre_phase = [phi];
ins_amplt = [A];
    
fTrue = cell(1,2);
numGroup = 1;
% fff = f1 + f2;

shapeTrue = cell(1,2);
tic;
[shapeInLoop,comp,errorRec,SL2Rec,iter,flag] = srcIterRegJC(x,N,numGroup,ins_amplt,ins_pre_phase,opt,fTrue,shapeTrue);
tiempo_RDBR = toc

opt.maxiter = 300;
opt.eps_error = 1e-6;
opt.show = 0;
opt.iterStyle = 'GS';
opt.shapeMethod = 2;
opt.eps_diff = 1e-6;
opt.ampErrBandWidth = 20;
opt.numSweep = 10;
    
switch opt.shapeMethod
    case 1
        opt.para.Ls=1000;
        opt.para.bandWidth = 10;
        opt.para.diffeoMethod = 'nufft';
    case 2
        opt.para.nknots = 10;
        opt.para.knotremoval_factor= 1.0001;
        opt.para.order = 3;
        opt.para.Ls = 1000;
end
tic;
inst_freq = zeros(1,N);
inst_freq(1,1:end-1) = N*diff(phi)*0.5/pi; inst_freq(1,end) = inst_freq(1,end-1);

[shape,component,Hcoef,flag,idx,iter,iterDR] = DeCom_MMD(x,t,1,[A],inst_freq,[phi],opt);
tiempo_MMD = toc
%% Extension 
Np = round(0.1*N);

s_ext = extendSig(x,phi,3,Np,'fw-bw');

fprintf('Extension Completed \n')

%% STFT of the extended signal

Next = length(s_ext);
fmax = 0.1;
te = 0:1/fs:(Next-1)/fs;
fe = 0:fs/Next:fs*fmax-fs/Next;
sigma = compute_sigma(s_ext,1);
be = 12;

[STFTe, sFe] = STFT_Gauss(s_ext,Next,sigma,fmax);
RTF = STFTe;
c1e = ridge_ext(RTF,0,0,10,10);
[A_ext,phi_ext] = extract_harmonics(RTF,sFe,c1e,be,be,1);
figure(2)
subplot(311)
plot(te,s_ext)
subplot(3,1,[2,3])
spectro(RTF,te,fe)
hold on
plot(te,fe(c1e),'r')
plot(te,fe(c1e-be),'b--')
plot(te,fe(c1e+be),'b--')
hold off

%% tvWSE
r_opt = D1;
Cext = makeC(A_ext,phi_ext,r_opt);

ve = ((Cext'*Cext)\Cext')*s_ext;

sen = s_ext./(ve(1)*A_ext)';

Cen = makeC(ones(1,Next),phi_ext,r_opt);

ven = ((Cen'*Cen)\Cen')*sen;

r_max = floor(0.5*N/max(c1e));
fmax = 1.2*(r_opt/r_max);

[Fen, sFen] = STFT_Gauss(sen,Next,sigma,fmax);
c1e = ridge_ext(Fen,0,0,10,10);
nv = 0;
vnv = NNodes(nv,r_opt,Fen,sFen,c1e,be,fs,0.9);

vh = Init_tvWSE(ven,vnv,r_opt,1,N,Next,0,1);

[lb,ub] = create_bounds(vnv,vh,r_opt,N,1);

method = 'pchip';

fprintf('Computing tvWSE on extended signal using %s. Nro of coefs : %i \n',method,numel(vh))

tic;
[s_tvwse_n,v_ie,eflag_tvwse] = tvWSE(sen',ones(1,Next),phi_ext,r_opt,vnv,vh,method,lb,ub,1,1,1);
t_tvwse = toc;

se_tvwse = ve(1)*A_ext'.*s_tvwse_n';
s_tvwse = se_tvwse(Np+1:length(x)+Np);
figure(3)
subplot(211)
plot(t,x)
hold on
plot(t,s_samd)
plot(t,s_tvwse)
hold off
legend('x','SAMD','tvWSE')

%% 
ind = 8*fs+1:28*fs;
x1 = x(ind) - mean(x(ind));

s_wse = s_tvwse(ind) - mean(s_tvwse(ind));
    
s_s = s_samd(ind) - mean(s_samd(ind));

s_md = component{1}(ind) - mean(component{1}(ind));
    
r_wse = x1-s_wse';
r_s = x1-s_s;
r_md = x1-s_md;

r_wse = r_wse - mean(r_wse);
r_md = r_md - mean(r_md);

Lag = 2000;
ccorr_wse = xcorr(s_wse,r_wse,Lag,'unbiased')/xcorr(s_wse,r_wse,0,'unbiased');
ccorr_s = xcorr(s_s,r_s,Lag,'unbiased')/xcorr(s_s,r_s,0,'unbiased');
ccorr_md = xcorr(s_md,r_md,Lag,'unbiased')/xcorr(s_md,r_md,0,'unbiased');

acorr_wse = xcorr(r_wse,Lag,'unbiased')/xcorr(r_wse,0,'unbiased');
acorr_s = xcorr(r_s,Lag,'unbiased')/xcorr(r_s,0,'unbiased');
acorr_md = xcorr(r_md,Lag,'unbiased')/xcorr(r_md,0,'unbiased');

R_wse = corrcoef(s_wse,r_wse);
R_wse = R_wse(1,2);

R_s = corrcoef(s_s,r_s);
R_s = R_s(1,2);
        
S_wse = abs((fft(acorr_wse)));
S_wse = S_wse(1:round(end/2));

S_s = abs((fft(acorr_s)));
S_s = S_s(1:round(end/2));

Sn_wse = S_wse/sum(S_wse);
SE_wse = -sum(Sn_wse.*log(Sn_wse));

Sn_s = S_s/sum(S_s);
SE_s = -sum(Sn_s.*log(Sn_s));

S_md = abs((fft(acorr_md)));
S_md = S_md(1:round(end/2));

%% Signal Plot 
fntsz = 12;
xi = x(ind) - mean(x(ind));
Ni = length(xi);
del = 1.96/sqrt(Ni);
figure(1);
set(gcf,'Position',[672 204 1245 758])
subplot(4,5,1:3);
plot(t(ind)-t(ind(1)),xi,'k')
text(0.05,900,'Recording 44; T6-O2 channel','FontSize',11)
ylim([-600 750])
xlim([t(1) t(length(ind))])
set(gca,'FontSize',fntsz)

col = [1 0 0];
subplot(4,5,6:8);
plot(t(ind)-t(ind(1)),xi,'k')
hold on
plot(t(ind)-t(ind(1)),s_tvwse(ind),'Color',col,'LineWidth',1.5);
hold off
text(0.05,580,'Ours','FontSize',12)
ylim([-500 500])
xlim([t(1) t(length(ind))])
set(gca,'FontSize',fntsz)

subplot(4,5,11:13);
plot(t(ind)-t(ind(1)),xi,'k')
hold on
plot(t(ind)-t(ind(1)),s_samd(ind),'Color',col,'LineWidth',1.5);
hold off
text(0.05,580,'SAMD','FontSize',12)
ylim([-500 500])
xlim([t(1) t(length(ind))])
set(gca,'FontSize',fntsz)

subplot(4,5,16:18);
plot(t(ind)-t(ind(1)),xi,'k')
hold on
plot(t(ind)-t(ind(1)),component{1}(ind),'Color',col,'LineWidth',1.5);
text(0.05,580,'MMD','FontSize',12)
hold off
ylim([-500 500])
xlim([t(1) t(length(ind))])
xlabel('Time [s]','FontSize',16)
set(gca,'FontSize',fntsz)

subplot(4,5,9:10)
plot(acorr_wse(Lag:end))
hold on
yline(del)
yline(-del)
xlim([0 Lag])
text(20,0.85,'Autocorr Res. Ours','FontSize',14)
hold off
set(gca,'FontSize',fntsz)

subplot(4,5,14:15)
plot(acorr_s(Lag:end))
hold on
yline(del)
yline(-del)
xlim([0 Lag])
text(20,0.85,'Autocorr Res. SAMD','FontSize',14)
hold off
set(gca,'FontSize',fntsz)

subplot(4,5,19:20)
plot(acorr_md(Lag:end))
hold on
yline(del)
yline(-del)
xlim([0 Lag])
xlabel('Lag [s]','FontSize',14)
text(20,0.85,'Autocorr Res. MMD','FontSize',14)
hold off
set(gca,'FontSize',fntsz)