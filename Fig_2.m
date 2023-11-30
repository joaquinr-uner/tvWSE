% This Script generates Fig. 2 from the paper "Fully Adaptive Time-Varying
% Wave-Shape Model: Applications in Biomedical Signal Processing"

addpath(genpath('auxiliary-functions'))
addpath(genpath('time-frequency-analysis'))

mInterp = 'pchip';
ext = 1;
demod = 0;
Crit = {'Wang'};
Cparams = struct('c',[6,8,10,12]);
%% Synthethic Signal

fs = 2000;
N = 2000;
Np = 0.1*N;
Next = N + 2*Np;
T = N/fs;
t = 0:1/fs:T-1/fs;
f = 0:fs/N:fs/2-fs/N;

D = 3;
alpha = zeros(D,N);

A = 0.1*sqrt(t+1);

phi = 40*t+5/(2*pi)*sin(2*pi*t);

alpha(1,:) = ones(1,N);
alpha(2,:) = 0.5 + 0.25*cos(2*pi*3*t);
alpha(3,:) = 0.3 + 0.25*cos(2*pi*4*t);

e = [1 2.005 2.995];

sigma = 1e-4;
fmax = 0.5;
b = round(3/pi*sqrt(sigma/2)*N);
outn = 1;

s_c = alpha(1,:).*cos(2*pi*phi);
for l=2:D
    s_c = s_c + alpha(l,:).*(cos(2*pi*e(l)*phi));
end
s_c = A.*s_c;

s_c = s_c - mean(s_c);

s = s_c;
s = s - mean(s);

%% Fundamental Component Estimation

[F,sF] = STFT_Gauss(s,N,sigma,fmax);
c = ridge_ext(F,0.1,0.1,10,10);

r_max = floor(0.5*N/max(c));

[A_est,phi_est] = extract_harmonics(F,sF,c,b,b,1);

D = order_opt(s',r_max,A_est,phi_est,Crit,Cparams,F);

if D<2
    D = 2;
end
C = makeC(A_est,phi_est,D);
v = ((C'*C)\C')*s';
sn = s./(v(1)*A_est);

[Fn,sFn] = STFT_Gauss(sn,N,sigma,fmax);

nv = 0;
vnv = NNodes(nv,D,Fn,sFn,c,b,fs,0.7);
fl = diff(phi_est)*fs;

cycls = 3;
s_ext = extendSig(s,phi_est,cycls,Np,'fw-bw');


[Fext,sFext] = STFT_Gauss(s_ext,Next,sigma,fmax);
cext = ridge_ext(Fext,0.1,0.1,10,10);
be = round(3/pi*sqrt(sigma/2)*Next);
[A_ext,phi_ext] = extract_harmonics(Fext,sFext,cext,be,be,1);

C = construct_dct(A_ext,phi_ext,D);
ve = ((C'*C)\C')*s_ext;

%% LR
se_lr = C*ve;
se_lr = se_lr';
tic;

%% SAMD
[se_samd,~,v_e] = my_SAMD(s_ext',A_ext,phi_ext*2*pi,D,1);
t_samd = toc;

B1e = ve(1)*A_ext;
se_norm = s_ext'./B1e;
se_norm = se_norm - mean(se_norm);

Cn = construct_dct(ones(1,Next),phi_ext,D);
vn = ((Cn'*Cn)\Cn')*se_norm';

[Fen,sFen] = STFT_Gauss(se_norm,Next,sigma,fmax);

%% Algorithm Initialization
%vn = zeros(1,length(vn));
vh = Init_tvWSE(vn,vnv,D,1,N,Next,0,outn);
[lb,ub] = create_bounds(vnv,vh,D,Next,outn);

options = optimoptions('lsqcurvefit','Algorithm','levenberg-marquardt',...
    'MaxFunctionEvaluations',400*length(vh),'Display','off');
tic;

%% tvWSE
[se_tvwse, v_ie] = tvWSE(se_norm,ones(1,Next),phi_ext,D,vnv,vh,mInterp,lb,ub,1,ext,outn);

se_tvwse = B1e.*se_tvwse;

t_tvwse = toc;
%% Posprocessing
[ti,ql,gamh,eh] = parse_coefs(Next,v_ie,D,vnv,0,1);
[a,~,Ale,alp] = compute_hafs(ti,ql,gamh,'pchip',1,B1e);
s_lr = se_lr(Np+1:Next-Np);
s_samd = se_samd(Np+1:Next-Np);
s_tvwse = se_tvwse(Np+1:Next-Np);
Al = Ale(:,Np+1:Next-Np);

alp2 = alp{1};
alp3 = alp{2};

t1 = ti{1};
t2 = ti{2};

q1 = ql{1};
q2 = ql{2};

%% Graphics

vec_c = [0 0 0; 0.494 0.184 0.556; 0.929 0.694 0.125; 0.85 0.325 0.098];
figure(1)
set(gcf,'Position',[563, 260, 619, 532])
subplot(3,3,1)
plot(t(0.1*N-1:0.2*N),s_c(0.1*N-1:0.2*N),'LineWidth',1.75,'Color',vec_c(1,:))
hold on
plot(t(0.1*N-1:0.2*N),s_lr(0.1*N-1:0.2*N),'Color',vec_c(2,:))
plot(t(0.1*N-1:0.2*N),s_samd(0.1*N-1:0.2*N),'Color',vec_c(3,:))
plot(t(0.1*N-1:0.2*N),s_tvwse(0.1*N-1:0.2*N),'Color',vec_c(4,:))
xlim([0.1 0.2])
ylim([-0.5 0.5])
xlabel('Time [s]')
ylabel('A.U.')
hold off
subplot(3,3,2)
plot(t(0.5*N-1:0.6*N),s_c(0.5*N-1:0.6*N),'LineWidth',1.75,'Color',vec_c(1,:))
hold on
plot(t(0.5*N-1:0.6*N),s_lr(0.5*N-1:0.6*N),'Color',vec_c(2,:))
plot(t(0.5*N-1:0.6*N),s_samd(0.5*N-1:0.6*N),'Color',vec_c(3,:))
plot(t(0.5*N-1:0.6*N),s_tvwse(0.5*N-1:0.6*N),'Color',vec_c(4,:))
xlim([0.5 0.6])
ylim([-0.5 0.5])
xlabel('Time [s]')
leg1=legend('s','LR','SAMD','Our method','Location','northoutside','Orientation','horizontal');
set(leg1,'Box','off')
set(leg1,'Position',[0.3 0.9592 0.4029 0.0319])
hold off
set(gca,'Position',[0.4108 0.7093 0.2093 0.2157])
subplot(3,3,3)
plot(t(0.8*N-1:0.9*N),s_c(0.8*N-1:0.9*N),'LineWidth',1.75,'Color',vec_c(1,:))
hold on
plot(t(0.8*N-1:0.9*N),s_lr(0.8*N-1:0.9*N),'Color',vec_c(2,:))
plot(t(0.8*N-1:0.9*N),s_samd(0.8*N-1:0.9*N),'Color',vec_c(3,:))
plot(t(0.8*N-1:0.9*N),s_tvwse(0.8*N-1:0.9*N),'Color',vec_c(4,:))
xlim([0.8 0.9])
ylim([-0.5 0.5])
xlabel('Time [s]')
hold off
subplot(3,3,[4,5,6])
plot(t,s_c,'-k')
hold on
plot(t,s_tvwse','-r')
xlabel('Time [s]')
ylabel('A.U.')
ylim([-0.3 0.4])
leg2=legend('s(t)','Our method','Orientation','horizontal');
set(leg2,'Box','off')
hold off

ax10=subplot(3,3,7);
plot(t,alpha(2,:),'k')
hold on
plot(t,a(2)*Al(1,:),'r--')
ti2 = ceil(t1)-Np;
ti2 = ti2(2:end-1);
stem(t(round(ti2)),a(2)*alp2(2:end-1),'mo')
xlabel('Time [s]')
ylabel('A.U.')
set(ax10,'Position',[0.13 0.13 0.35 0.2])
ylim([0 1.1])
leg10=legend({'$\alpha_2(t)$','$\hat{\alpha}_2(t)$'},'Interpreter','latex','Orientation','horizontal');
set(leg10,'Box','off')
hold off

ax11=subplot(3,3,9);
plot(t,alpha(3,:),'k')
hold on
plot(t,a(3)*Al(2,:),'r--')
ti3 = ceil(t2)-Np;
ti3 = ti3(2:end-1);
stem(t(round(ti3)),a(3)*alp3(2:end-1),'mo')
xlabel('Time [s]')
ylim([0 0.85])
set(ax11,'Position',[0.56 0.13 0.35 0.2])
leg11=legend({'$\alpha_3(t)$','$\hat{\alpha}_3(t)$'},'Interpreter','latex','Orientation','horizontal');
set(leg11,'Box','off')
hold off