s = zpk('s'); 
G = (s-1)/(s+1);
W1 = 0.1*(s+100)/(100*s+1); 
W2 = 0.1; 
W3 = [];
P = augw(G,W1,W2,W3);

opt.AutoScale = 'off';
opt.Regularize = 'off';
[K,CL,GAM] = hinfsyn(P,1,1,0.187); 
[K2,CL2,GAM2] = h2syn(P,1,1,opt);

GAM
GAM2

L = G*K; 
S = inv(1+L); 
T = 1-S; 

KM = K;
K2M = K2;

sigma(S,'k',GAM/W1,'k-.',T,'r',GAM*G/W2,'r-.')
legend('S = 1/(1+L)','GAM/W1','T=L/(1+L)','GAM*G/W2','Location','Northwest')
shg

%%
s = zpk('s'); 
G = (s-1)/(s+1);
W1 = 0.1*(s+100)/(100*s+1); 
W2 = 0.1; 
W3 = [];
P = augw(G,W1,W2,W3);

[K,CL,GAM] = hinfgsyn(P,1,1,0.187); 
[K2,CL2,GAM2] = h2gsyn(P,1,1);

GAM
GAM2

L = G*K; 
S = inv(1+L); 
T = 1-S; 

KD = K;
K2D = K2;

sigma(S,'k',GAM/W1,'k-.',T,'r',GAM*G/W2,'r-.')
legend('S = 1/(1+L)','GAM/W1','T=L/(1+L)','GAM*G/W2','Location','Northwest')