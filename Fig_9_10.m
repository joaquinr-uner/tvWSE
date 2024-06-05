% This Script generates Figs. 9 and 10 from the paper "Fully Adaptive Time-Varying
% Wave-Shape Model: Applications in Biomedical Signal Processing"

addpath(genpath('time-frequency-analysis'))
addpath(genpath('auxiliary-functions'))

load('cu03m.mat')
st=1.12e5;
ed=1.18e5-1;

fs = 250;

s = val(st:ed);
s = s - mean(s);

N = length(s);

t = 0:1/fs:(N-1)/fs;

fmax = 0.1;

f = 0:fs/N:fs*fmax-fs/N;
%%
sigma = 5e-5;

b = round(3/pi*sqrt(sigma/2)*N);
[F, sF ] = STFT_Gauss(s,N,sigma,fmax);

c = ridge_ext(F,0.1,0.1,10,10);

c = ridge_correct(c,F,b);


figure(1)
subplot(311)
plot(t,s)
subplot(3,1,[2,3])
spectro(F,t,f)
hold on
plot(t,f(c),'r')

[A,phi ] = extract_harmonics(F,sF,c,b,b,1);
%% SAMD and tvWSE
r_opt = 6;
C = makeC(A,phi,r_opt);

v = (((C')*C)\C')*s';

s_lr = C*v;

[s_samd,~,v_samd] = my_SAMD(s,A,phi*(2*pi),r_opt,1);
fprintf('SAMD Completed \n')

figure(1)
subplot(311)
hold on
plot(t,s_samd)
hold off

Np = round(0.1*N);
s_ext = extendSig(s,phi,3,Np,'fw-bw');

Next = N + 2*Np;
fe = 0:fs/Next:fs*fmax-fs/Next;
[Fext, sFext] = STFT_Gauss(s_ext,Next,sigma,fmax);

ce = ridge_ext(Fext,0,0,10,10);
[A_ext,phi_ext] = extract_harmonics(Fext,sFext,ce,b,b,1);

Ce = makeC(A_ext,phi_ext,r_opt);

ve = (((Ce')*Ce)\Ce')*s_ext;

sen = s_ext./(ve(1)*A_ext');

Cn = makeC(ones(1,Next),phi_ext,r_opt);

vn = ((Cn'*Cn)\Cn')*sen;

r_max = floor(0.5*Next/max(ce));
fmax = 1.2*(r_opt/r_max);

[Fen, sFen] = STFT_Gauss(sen,Next,sigma,fmax);
nv = 0;
vnv = NNodes(nv,r_opt,Fen,sFen,ce,b,fs,0.9);

vh = Init_tvWSE(vn,vnv,r_opt,1,N,Next,0,1);

[lb,ub] = create_bounds(vnv,vh,r_opt,N,1);

method = 'pchip';

fprintf('Computing tvWSE on extended signal using %s. Nro of coefs : %i \n',method,numel(vh))

tic;
[s_tvwse_n,v_ie,eflag_tvwse] = tvWSE(sen',ones(1,Next),phi_ext,r_opt,vnv,vh,method,lb,ub,1,1,1);
t_tvwse = toc

se_tvwse = v(1)*A_ext.*s_tvwse_n;

s_tvwse = se_tvwse(0.1*N+1:N+0.1*N);

figure(3)
subplot(211)
plot(t,s)
hold on
plot(t,s_lr)
plot(t,s_samd)
plot(t,s_tvwse)
hold off

%% MMD

opt.maxiter = 10000;
opt.eps_error = 1e-12;
opt.show = 0;
opt.nknots = 20;
opt.knotremoval_factor= 1.01;
opt.order = 3;
opt.eps_diff = opt.eps_error;
    
ins_pre_phase = [phi]*2*pi;
ins_amplt = [A];
    
fTrue = cell(1,2);
numGroup = 1;
% fff = f1 + f2;

shapeTrue = cell(1,2);
%         shapeTrue{1} = @(x) sh1(x);
%         shapeTrue{2} = @(x) sh2(x);
tic;
[shapeInLoop,comp,errorRec,SL2Rec,iter,flag] = srcIterRegJC(s,N,numGroup,ins_amplt,ins_pre_phase,opt,fTrue,shapeTrue);
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
inst_freq(1,1:end-1) = N*diff(phi); inst_freq(1,end) = inst_freq(1,end-1);

[shape,component,Hcoef,flag,idx,iter,iterDR] = DeCom_MMD(s,t,1,[A],inst_freq,[phi],opt);
tiempo_MMD = toc

%% Compute HAFs
Next = 1.2*N;
[ti, alp, gamh, eh] = parse_coefs(Next,v_ie,r_opt,vnv,0);
[ti,alp] = remove_outn(ti,alp,1);
[a, b, Q, q] = compute_hafs(ti,alp,gamh,'pchip',1,A);

%% Fig. 9
t = 0:1/fs:length(s)/fs-1/fs;
ind1 = 2*fs+1:6*fs;

ind2 = 11*fs+1:15*fs;
ind3 = 20*fs+1:24*fs;
col = 'k';
fonts = 13;
figure(3)
subplot(5,3,[1,2,3])
plot(t,s,'Color',col)
hold on
ylim([-750 750])
yl = ylim;
offs = 50;
patch([2,6,6,2],[yl(2)-offs yl(2)-offs yl(1)+offs yl(1)+offs],'red','FaceColor','none','EdgeColor','red','LineStyle','--')
patch([11,15,15,11],[yl(2)-offs yl(2)-offs yl(1)+offs yl(1)+offs],'red','FaceColor','none','EdgeColor','red','LineStyle','--')
patch([20,24,24,20],[yl(2)-offs yl(2)-offs yl(1)+offs yl(1)+offs],'red','FaceColor','none','EdgeColor','red','LineStyle','--')
title('ECG Signal with ventricular fibrillation event')
set(gca,'FontSize',fonts)
xlim([t(1) t(end)])
subplot(5,3,4)
plot(t(ind1),s(ind1),'Color',col)
ylabel('Signal')
set(gca,'FontSize',fonts)
set(gca,'YTick',[-200 600])
%title('Signal')
subplot(5,3,7)
plot(t(ind1),s_tvwse(ind1),'Color',col)
ylabel('Ours')
set(gca,'FontSize',fonts)
set(gca,'YTick',[-200 400])
%title('Ours')
subplot(5,3,10)
plot(t(ind1),s_samd(ind1),'Color',col)
ylabel('SAMD')
set(gca,'FontSize',fonts)
set(gca,'YTick',[-100 200])
%title('SAMD')
subplot(5,3,13)
plot(t(ind1),component{1}(ind1),'Color',col)
ylabel('MMD')
set(gca,'FontSize',fonts)
set(gca,'YTick',[-100 400])
%title('MMD')
%subplot(6,3,16)
%plot(t(ind1),s_lr(ind1))
%set(gca,'FontSize',fonts)
%set(gca,'YTick',[-100 200])
%ylabel('LR')
xlabel('Time [s]')
%title('LR')
subplot(5,3,5)
plot(t(ind2),s(ind2),'Color',col)
set(gca,'FontSize',fonts)
set(gca,'YTick',[-200 400])
%title('Signal')
subplot(5,3,8)
plot(t(ind2),s_tvwse(ind2),'Color',col)
set(gca,'FontSize',fonts)
set(gca,'YTick',[-200 400])
%title('Ours')
subplot(5,3,11)
plot(t(ind2),s_samd(ind2),'Color',col)
set(gca,'FontSize',fonts)
set(gca,'YTick',[-200 200])
%title('SAMD')
subplot(5,3,14)
plot(t(ind2),component{1}(ind2),'Color',col)
set(gca,'FontSize',fonts)
set(gca,'YTick',[-200 400])
%title('MMD')
%subplot(6,3,17)
%plot(t(ind2),s_lr(ind2))
%set(gca,'FontSize',fonts)
%set(gca,'YTick',[-200 200])
xlabel('Time [s]')
%title('LR')
subplot(5,3,6)
plot(t(ind3),s(ind3),'Color',col)
set(gca,'FontSize',fonts)
set(gca,'YTick',[-500 500])
set(gca,'YLim',[-663 626])
%title('Signal')
subplot(5,3,9)
plot(t(ind3),s_tvwse(ind3),'Color',col)
set(gca,'FontSize',fonts)
set(gca,'YTick',[-500 500])
%title('Ours')
subplot(5,3,12)
%title('SAMD')
plot(t(ind3),s_samd(ind3),'Color',col)
set(gca,'FontSize',fonts)
subplot(5,3,15)
plot(t(ind3),component{1}(ind3),'Color',col)
set(gca,'YTick',[-1000 1000])
set(gca,'FontSize',fonts)
%title('MMD')
%subplot(6,3,18)
%plot(t(ind3),s_lr(ind3))
%set(gca,'FontSize',fonts)
%set(gca,'YTick',[-500 500])
xlabel('Time [s]')
%title('LR')

%% Fig. 10
figure(5)
set(gcf,'Position',[680 488 1241 474])
subplot(3,1,[1,2])
plot(t,s);
hold on
scatter(t,s,5,a(2)*Q(1,:),'filled','Marker','o')
set(gca,'FontSize',14)
ylim([-750 750])
xlim([t(1) t(end)])
ylabel('U.A.')
%title('ECG with ventricular fibrillation event')
title('ECG con evento de fibrilación ventricular')
%title('ECG Signal with ventricular fibrillation event')
colorbar 
subplot(3,1,3)
fill([t fliplr(t)],[a(2)*Q(1,:)-0.05 fliplr(a(2)*Q(1,:)+0.05)],[a(2)*Q(1,:) fliplr(a(2)*Q(1,:))])
%title('Second Harmonic Amplitude Function')
title('Segunda Función de Amplitud Armónica')
colorbar
set(gca,'FontSize',14)
xlim([t(1) t(end)])
xlabel('Tiempo [s]')
ylabel('U.A.')
