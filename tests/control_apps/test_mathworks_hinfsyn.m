close all

s = zpk('s');
G = (s-1)/(s+1);
W1 = 0.1*(s+100)/(100*s+1); 
W2 = 0.1; 
W3 = [];
P = augw(G,W1,W2,W3);

% g = .42;
% g = .2;
g = 0.1847;
% g = 0.18278;
[KM,CLM,gammaM] = hinfsyn(P,1,1,g);
[K,CL,gamma] = hinfgsyn(P,1,1,g);

gammaM
gamma

% eig(CLM.a)
% eig(CL.a)

% figure
% sigma(CLM,CL,ss(g)); grid
% legend('H\infty MATLAB','H\infty defsub','\gamma');

% figure
% bode(KM,K); grid
% legend('H\infty MATLAB','H\infty defsub');

figure
sigmaplot(CLM,CL,ss(g)); grid
legend('H\infty MATLAB','H\infty defsub','\gamma');
