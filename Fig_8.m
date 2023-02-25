addpath(genpath('/home/sentey/Dropbox/Github'))

%% Data Loading
x = load('2869156_Resp(32).txt');
ecg = load('2869156_II(256).txt');
pleth = load('2869156_Pleth(64).txt');

x = x(:)';
fs = 32;
x = x(1:fs*60*10);
N = length(x);
t = 0:1/fs:N/fs-1/fs;

fs_ecg = 256;
ecg = ecg(1:fs_ecg*60*10);
t_ecg = 0:1/fs_ecg:length(ecg)/fs_ecg-1/fs_ecg;

fs_pleth = 64;
pleth = pleth(1:fs_pleth*60*10)';
t_pleth = 0:1/fs_pleth:length(pleth)/fs_pleth-1/fs_pleth;

s = x - mean(x);
s = s';

%% LR and SAMD

%Parameters for STFT
bt = 1:N;
redun = 1;
gamma = 0;
sigma = 0.045;
fmax = 0.05;
ft = 1/redun:1/redun:round(fmax*N);

[STFT,SST1,SST2] = sstn_test_modL_new(s',gamma,sigma,ft,bt,redun);

RTF = SST2;
d = 2;
c1 = exridge(RTF,0,0,2);
aux = zeros(size(t));
for j = 1:N
    aux(j) = sum(RTF(c1(j)-d:c1(j)+d,j));
    RTF(c1(j)-d:c1(j)+d,j) = 0;
end
A1_est = abs(aux); 
phi1_est = phase(aux)/(2*pi);

c2 = exridge(RTF(201:350,:),0,0,2)+200;
aux = zeros(size(t));
for j = 1:N
    aux(j) = sum(RTF(c2(j)-d:c2(j)+d,j));
end
A2_est = abs(aux); 
phi2_est = phase(aux)/(2*pi);

D1 = 2; 
D2 = 5; 

tic;
C = [makeC(A1_est,phi1_est,D1) makeC(A2_est,phi2_est,D2)];

C1_aux = makeC(ones(1,1000),-pi:2*pi/1000:pi-2*pi/1000,D1);
C2_aux = makeC(ones(1,1000),-pi:2*pi/1000:pi-2*pi/1000,D2);

coef_est = ((C'*C)\C')*s;
f_est_LR = C*coef_est;
f1_est_LR = C(:,1:2*D1)*coef_est(1:2*D1);
f2_est_LR = C(:,2*D1+1:end)*coef_est(2*D1+1:end);
tiempo_LR = toc

WSF1 = C1_aux*coef_est(1:2*D1);
WSF2 = C2_aux*coef_est(2*D1+1:end);

order_Phi = 2;
warning('off')
tic;
[f_est_NLR,f1_est_NLR,f2_est_NLR,b_e] = two_wsf_nonlinear(s',A1_est,phi1_est*2*pi,D1,A2_est,phi2_est*2*pi,D2,order_Phi);
tiempo_NLR = toc
warning('on')

%% RDBR and MMD

opt.maxiter = 10000;
opt.eps_error = 1e-12;
opt.show = 0;
opt.nknots = 20;
opt.knotremoval_factor= 1.01;
opt.order = 3;
opt.eps_diff = opt.eps_error;
    
ins_pre_phase = [phi1_est;phi2_est];
ins_amplt = [A1_est;A2_est];
    
fTrue = cell(1,2);
numGroup = 2;

shapeTrue = cell(1,2);
         shapeTrue{1} = @(x) sh1(x);
         shapeTrue{2} = @(x) sh2(x);%
tic;
[shapeInLoop,comp,errorRec,SL2Rec,~,~] = srcIterRegJC(s',N,numGroup,ins_amplt,ins_pre_phase,opt,fTrue,shapeTrue);
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
 inst_freq = zeros(2,N);
 inst_freq(1,1:end-1) = N*diff(phi1_est); inst_freq(1,end) = inst_freq(1,end-1);
 inst_freq(2,1:end-1) = N*diff(phi2_est); inst_freq(2,end) = inst_freq(2,end-1);
 [shape,component,Hcoef,flag,idx,iter,iterDR] = DeCom_MMD(s',t,2,[A1_est;A2_est],inst_freq,[phi1_est;phi2_est],opt);
 tiempo_MMD = toc

%% Signal Extension
s = x;
s = s - mean(s);

phi = phi1_est;
A = A1_est;
N = length(s);
Np = round(0.1*N);

s_ext = extendSig(s,phi1_est/(2*pi),3,Np,'fw-bw');
fprintf('Extension Completed \n')
fmax = 0.05;
Next = N + 2*Np;
te = 0:1/fs:(Next-1)/fs;
fe = 0:fs/Next:fs*fmax-fs/Next;
%% tvWSE 1st Comp.
sigma = 1e-6;
be = round(3/pi*sqrt(sigma/2)*Next);

[Fext, sFext ] = STFT_Gauss(s_ext,Next,sigma,fmax);
Uext = istct_fast(Fext,fe,0.3);
Wext = Fext.*Uext;

ce = ridge_ext(Fext,0.1,0.1,10,10);
ce = ridge_correct(ce,Fext,b);

figure(2)
subplot(311)
plot(te,s_ext)
subplot(3,1,[2,3])
spectro(Fext,te,fe)
hold on
plot(te,fe(ce),'r')
plot(te,fe(ce-be),'b--')
plot(te,fe(ce+be),'b--')
hold off

[A_ext, phi_ext] = extract_harmonics(Fext,sFext,ce,be,be,1);

Cext = construct_dct(A_ext,phi_ext,r_opt);

ve = ((Cext'*Cext)\Cext')*s_ext;

sen = s_ext./(ve(1)*A_ext');

Cen = construct_dct(ones(1,Next),phi_ext,r_opt);

ven = ((Cen'*Cen)\Cen')*sen;

fmax = 0.05;

[Fen, sFen] = STFT_Gauss(sen,Next,sigma,fmax);
nv = 0;
vnv = NNodes(nv,D1,Fen,sFen,ce,be,fs,0.9);

vh = Init_tvWSE(ven,vnv,D1,1,N,Next,0,1);

[lb,ub] = create_bounds(vnv,vh,D1,N,1);

method = 'pchip';

fprintf('Computing tvWSE on extended signal using %s. Nro of coefs : %i \n',method,numel(vh))

tic;
[s_tvwse_n,v_ie,eflag_tvwse] = tvWSE(sen',ones(1,Next),phi_ext,D1,vnv,vh,method,lb,ub,1,1,1);
t_tvwse = toc;

se_tvwse = ve(1)*A_ext.*s_tvwse_n;

s_tvwse = se_tvwse(Np+1:N+Np);
%% tvWSE 2nd Comp.

s2i = s - s_tvwse;

sigma2 = compute_sigma(s2i,1);
b2 = round(3/pi*sqrt(sigma2/2)*N);
fmax2 = 0.25;

f = 0:fs/N:fs*fmax2-fs/N;
c2 = ridge_ext(F(201:350,:),0.1,0.1,10,10)+200;

[A2,phi2] = extract_harmonics(F2,sF2,c2,b2,b2,1);

r2_opt = 5;
A1 = A;
phi1 = phi;
r1_opt = r_opt;
C2 = makeC(A2,phi2,r2_opt);

v2 = ((C2'*C2)\C2')*s2l ;

s2_lr = C2*v2;

[~,modes_samd,v_samd] = my_SAMD(s,[A1_est;A2_est],[phi1_est;phi2_est],[D1,D2],1);
fprintf('SAMD Completed \n')

s2_ext = extendSig(s2i,phi2,3,0.1*N,'fw-bw');
fprintf('Extension Completed \n')
fmax = 0.05;
Next = N + 0.2*N;
te = 0:1/fs:(Next-1)/fs;
fe = 0:fs/Next:fs*fmax-fs/Next;

be = round(3/pi*sqrt(sigma/2)*Next);

[F2e, sF2e ] = STFT_Gauss(s2_ext,Next,sigma,fmax);

c2e = ridge_ext(F2e(280:420,:),0.1,0.1,10,10)+280;
%c2e = ridge_correct(c2e,F2e,b);


[A2e, phi2e] = extract_harmonics(F2e,sF2e,c2e,be,be,1);

C2e = makeC(A2e,phi2e,r2_opt);

v2e = ((C2e'*C2e)\C2e')*s2_ext;

s2en = s2_ext./(v2e(1)*A2e');

C2en = makeC(ones(1,Next),phi2e,r2_opt);

v2en = ((C2en'*C2en)\C2en')*s2en;

r_max = floor(0.5*N/max(c2));
fmax2 = 1.2*(r2_opt/r_max);

[F2en, sF2en] = STFT_Gauss(s2en,Next,sigma,fmax2);
nv = 0;
vnv2 = NNodes(nv,r2_opt,F2en,sF2en,c2e,be,fs,0.8);

vh2 = Init_tvWSE(v2en,vnv2,r2_opt,1,N,Next,0,1);

[lb2,ub2] = create_bounds(vnv2,vh2,r2_opt,N,1);

method = 'pchip';

fprintf('Computing tvWSE on extended signal using %s. Nro of coefs : %i \n',method,numel(vh2))

tic;
[s2_tvwse_n,v2_ie,eflag2_tvwse] = tvWSE(s2en',ones(1,Next),phi2e,r2_opt,vnv2,vh2,method,lb2,ub2,1,1,1);
t2_tvwse = toc

s2e_tvwse = v2e(1)*A2e.*s2_tvwse_n;

s2_tvwse = s2e_tvwse(Np+1:N+Np);


%% 
ecg = ecg - mean(ecg);
yl_e = [-200 350];
ind_e = 1:length(t_ecg);
figure(2);
h21= subplot(5,1,1);
plot(t_ecg(ind_e)-t_ecg(ind_e(1)),ecg(ind_e),'k')
text(0.05,1.2*yl_e(2),'Electrocardiogram')
ylim(yl_e)
xlim([t_ecg(1) t_ecg(length(ind_e))])

yl = [-200 200];
ecg2 = ecg/2;
h22 = subplot(5,1,2);
plot(t_ecg(ind_e)-t_ecg(ind_e(1)),ecg2(ind_e),'r--')
hold on
plot(t(ind)-t(ind(1)),s2_tvwse(ind),'g-','LineWidth',1.5);
hold off
text(0.05,1.2*yl(2),'tvWSE')
ylim(yl)
xlim([t(1) t(length(ind))])

h23 = subplot(5,1,3);
plot(t_ecg(ind_e)-t_ecg(ind_e(1)),ecg2(ind_e),'r--')
hold on
plot(t(ind)-t(ind(1)),f2_est_NLR(ind),'g-','LineWidth',1.5);
hold off
text(0.05,1.2*yl(2),'SAMD')
ylim(yl)
xlim([t(1) t(length(ind))])

linkaxes([h21 h22 h23],'x')
%%
    
sc_tvwse = s_tvwse + s2_tvwse;
sc_samd = f1_est_NLR + f2_est_NLR;

%%
figure(3);
yl = [-2000 2000];
ind = 1:N;
h31= subplot(5,1,1);
plot(t(ind)-t(ind(1)),s(ind),'k')
text(0.05,1.2*yl(2),'Impedance Pneumography')
ylim(yl)
xlim([t(1) t(length(ind))])


h32 = subplot(5,1,2);
plot(t(ind)-t(ind(1)),s(ind),'k')
hold on
plot(t(ind)-t(ind(1)),sc_tvwse(ind),'g-','LineWidth',1.5);
hold off
text(0.05,1.2*yl(2),'tvWSE')
ylim(yl)
xlim([t(1) t(length(ind))])

h33 = subplot(5,1,3);
plot(t(ind)-t(ind(1)),s(ind),'k')
hold on
plot(t(ind)-t(ind(1)),sc_samd(ind),'g-','LineWidth',1.5);
hold off
text(0.05,1.2*yl(2),'SAMD')
ylim(yl)
xlim([t(1) t(length(ind))])

linkaxes([h31 h32 h33],'x')
%%

tin = 185;
tend = 205;
ind = tin*fs:tend*fs;
ind_ecg = tin*fs_ecg:tend*fs_ecg;

t_m = linspace(min(phi2(ind)),max(phi2(ind)),N);
N = length(f2_est_NLR(ind));
T = length(WSF2);
P = (t_m(end)-t_m(1)+1);
P = floor(P);
m = interp1(phi2(ind),1:N,t_m,'spline');

y = interp1(1:N,s2_tvwse(ind)./A2(ind),m,'spline');
waves_tvwse = zeros(P,T); 
for i = 1:P
    waves_tvwse(i,:) = interp1(t_m,y,t_m(1)+i-1+[0:1/(T):1-1/(T)],'spline');
end

t_m = linspace(min(phi2_est(ind)),max(phi2_est(ind)),N);
P = ((t_m(end)-t_m(1)+1)/(2*pi));
P = floor(P);
m = interp1(phi2_est(ind),1:N,t_m,'spline');

y = interp1(1:N,f2_est_NLR(ind)./A2_est(ind),m,'spline');
waves_SAMD = zeros(P,T); 
for i = 1:P
    waves_SAMD(i,:) = interp1(0.5*t_m/pi,y,0.5*t_m(1)/pi+i-1+[0:1/(T):1-1/(T)],'spline');
end

y = interp1(1:N,component{2}(ind)./A2_est(ind),m,'spline');
waves_MMD = zeros(P,T);
for i = 1:P
    waves_MMD(i,:) = interp1(0.5*t_m/pi,y,0.5*t_m(1)/pi+i-1+[0:1/(T):1-1/(T)],'spline');
end
fnts = 11;
figure(3);
set(gcf,'Position',[2231 340 1176 575])
subplot(4,9,1:4);
plot(t(ind),x(ind),'k')
text(tin,3550,'IP recording','FontSize',12,'FontWeight','bold')
hold off
%ylim([-10 15])
yticks([1000:1000:3000])
xlim([t(ind(1)) t(ind(end))])
xticks([])
set(gca,'Position',[0.11 0.7673 0.3326 0.1577])
set(gca,'FontSize',fnts)

subplot(4,9,10:13)
plot(t(ind),s_tvwse(ind),'k')
text(tin,1200,'Ours','FontSize',12)
text(tin+5,1200,'Respiratory Component','FontSize',12,'FontWeight','bold')
hold off
% ylim([-6.5 6.5])
xlim([t(ind(1)) t(ind(end))])
xticks([])
set(gca,'YTick',[-500 500])
set(gca,'Position',[0.11 0.5482 0.3326 0.1577])
set(gca,'FontSize',fnts)

subplot(4,9,14:17)
plot(t_ecg(ind_ecg),ecg(ind_ecg),'r-.');hold on; plot(t(ind),s2_tvwse(ind),'k'); 
text(tin,490,'Ours','FontSize',12)
text(tin+5,490,'Cardiac Component','FontSize',12,'FontWeight','bold')
hold off
ylim([-250 400])
xlim([t(ind(1)) t(ind(end))])
xticks([])
set(gca,'YTick',[-200 400])
set(gca,'FontSize',fnts)

subplot(4,9,18)
for i=1:P
plot(linspace(-0.5,0.5,length(WSF2)),waves_tvwse(i,:)); hold on
end
text(-0.35,4,{'Cardiac','Cycles'},'FontWeight','bold')
hold off
set(gca,'FontSize',fnts)

subplot(4,9,19:22)
plot(t(ind),f1_est_NLR(ind),'k')
text(tin,1050,'SAMD','FontSize',12)
hold off
% ylim([-6.5 6.5])
xlim([t(ind(1)) t(ind(end))])
set(gca,'YTick',[-500 500])
xticks([])
set(gca,'Position',[0.11 0.3291 0.3326 0.1577])
set(gca,'FontSize',fnts)

subplot(4,9,23:26)
plot(t_ecg(ind_ecg),ecg(ind_ecg),'r-.');hold on; plot(t(ind),f2_est_NLR(ind),'k'); 
text(tin,490,'SAMD','FontSize',12)
hold off
ylim([-250 400])
xlim([t(ind(1)) t(ind(end))])
xticks([])
set(gca,'YTick',[-200 400])
set(gca,'FontSize',fnts)

subplot(4,9,27)
for i=1:P
plot(linspace(-0.5,0.5,length(WSF2)),waves_SAMD(i,:)); hold on
end
set(gca,'FontSize',fnts)

subplot(4,9,28:31)
plot(t(ind),component{1}(ind),'k');
text(tin,1150,'MMD','FontSize',12)
hold off
% ylim([-11 11])
xlim([t(ind(1)) t(ind(end))])
xlabel('Time [s]','FontSize',12)
set(gca,'XTick',[tin:5:tend])
set(gca,'YTick',[-500 500])
set(gca,'Position',[0.11 0.11 0.3326 0.1577])
set(gca,'FontSize',fnts)

subplot(4,9,32:35)
plot(t_ecg(ind_ecg),ecg(ind_ecg),'r-.');hold on; plot(t(ind),component{2}(ind),'k');
text(tin,700,'MMD','FontSize',12)
hold off
ylim([-250 600])
xlim([t(ind(1)) t(ind(end))])
xlabel('Time [s]','FontSize',14)
set(gca,'XTick',[tin:5:tend])
set(gca,'YTick',[-200 600])
set(gca,'FontSize',fnts)

subplot(4,9,36)
for i=1:P
plot(linspace(-0.5,0.5,length(WSF2)),waves_MMD(i,:)); hold on
end
ylim([-0.6 1])
set(gca,'FontSize',fnts)
% set(gca,'yticklabel',[])