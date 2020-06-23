%% continuous-time system
load hinfsynExData P
ncont = 1; 
nmeas = 2; 

% P.d = rand(5,4);

opts = hinfsynOptions('Display','on');

% gamma = 0.9396;
gamma = 0.9397;

% [KM,CLM,GAMM,infoM] = hinfsyn(P,nmeas,ncont,0.9397,opts);
[KM,CLM,GAMM,infoM] = hinfsyn(P,nmeas,ncont,[gamma,1.45],opts);
[K,CL,GAM,info] = hinfgsyn(P,nmeas,ncont,gamma);

GAMM
GAM

bounded_real_lemma(CL.a,CL.b,CL.c,CL.d,3,false)

% bode(CLM,CL)
% legend('H\infty MATLAB','H\infty defsub');

% step(CLM,CL)
% bode(KM,K)
% legend('H\infty MATLAB','H\infty defsub');

% sigma(CLM,CL,ss(gamma),'-*',logspace(1,2,2000)); grid % matlab nu e sub
% gamma defapt...
sigma(CLM,CL,ss(gamma)); grid
% legend('H\infty MATLAB','H\infty defsub','\gamma');
legend(['Hinf MATLAB, \gamma = ',num2str(GAMM)],['Hinf defsub, \gamma = ',num2str(GAM)],['\gamma = ',num2str(gamma)]);

%% discrete-time system
load hinfsynExData P
ncont = 1; 
nmeas = 2; 

Pd = c2d(P,1e-2,'tustin');
% step(P,Pd)

P = Pd;
gamma = 0.94;

opts = hinfsynOptions('Display','on');
[KM,CLM,GAMM,INFO_M] = hinfsyn(P,2,1,gamma);
% pole(CLM)

[K,CL,GAM,INFO] = hinfgsyn(P,2,1,gamma);
% pole(CL)

GAMM
GAM

% step(CLM,CL,'*')
% legend('H\infty MATLAB','H\infty defsub');

% bode(KM,K)
% legend('H\infty MATLAB','H\infty defsub');

sigma(CLM,CL,ss(gamma),logspace(-2,3,1000)); grid
legend('H\infty MATLAB','H\infty defsub','\gamma');
