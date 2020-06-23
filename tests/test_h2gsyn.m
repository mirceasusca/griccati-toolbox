rng(1,'twister');
P = rss(2,4,5)';
P.A = P.A + 20*eye(size(P.A)); % cannot stabilize system
% P.D(1:3,1:3) = 0*P.D(1:3,1:3);
pole(P)

opt = rctoptions.h2syn;
opt.AutoScale = 'off';
opt.Regularize = 'off';
[KM,CLM,GAMM] = h2syn(P,2,1,opt);
% pole(CL)

[K,CL,GAM,INFO] = h2gsyn(P,2,1);
% pole(CLM)

GAMM
GAM 

% figure
% sigma(CLM,CL,ss(GAM))
% legend('H2 MATLAB','H2 defsub');

% step(CLM,CL)
sigmaplot(CLM,CL,ss(GAMM),ss(GAM)), grid
legend(['H2 MATLAB, \gamma = ',num2str(GAMM)],['H2 defsub, \gamma = ',num2str(GAM)],'GAMM','GAM');

%%
rng(14,'twister');
P = drss(2,5,4);
P.Ts = 0.01;
% P.A = P.A + eye(size(P.A)); % cannot stabilize system
% P.D(1:3,1:3) = zeros(3);
P.D(1:3,1:3) = rand(3);
pole(P)

opt = rctoptions.h2syn;
opt.AutoScale = 'off';
opt.Regularize = 'off';
[KM,CLM,GAMM] = h2syn(P,2,1,opt);
% pole(CLM)

[K,CL,GAM,INFO] = h2gsyn(P,2,1);
% pole(CL)

GAMM
GAM

% step(CLM,CL); grid;
% bode(KM,K), grid
sigmaplot(CLM,CL,ss(GAMM),ss(GAM)), grid
% sigma(CLM,CL,ss(GAMM),ss(GAM)), grid
legend(['H2 MATLAB, \gamma = ',num2str(GAMM)],['H2 defsub, \gamma = ',num2str(GAM)],'GAMM','GAM');