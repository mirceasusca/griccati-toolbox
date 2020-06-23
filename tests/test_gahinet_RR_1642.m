%% Example 1

A = 1;
B1 = [1,0];
C1 = [1;0];
B2 = 1;
C2 = 1;
D12 = [0;1];
D21 = D12';
D11 = zeros(2);
D22 = 0;

B = [B1,B2];
C = [C1;C2];
D = [D11,D12;D21,D22];

P = ss(A,B,C,D);
nmeas = 1;
ncont = 1;

g = 2.74;
[KM,CLM,GAMM] = hinfsyn(P,nmeas,ncont,g);
[K,CL,GAM] = hinfgsyn(P,nmeas,ncont,g);

GAMM
GAM

sigma(CLM,CL,ss(g));
legend('matlab','defsub','\gamma');

%% Example 2
A = [1,-1,0; 1,1,-1; 0,1,-2];
B1 = [1,2,0;0,-1,0;1,1,0];
B2 = [1;0;1];
C1 = [0,0,0;1,1,0;-1,0,1];
C2 = [0,-1,1];
D12 = [1;0;0];
D21 = [0,0,1];
D11 = zeros(3);
D22 = 1;

B = [B1,B2];
C = [C1;C2];
D = [D11,D12;D21,D22];

P = ss(A,B,C,D);
nmeas = 1;
ncont = 1;

g = 30;
% g =21.527876;
[KM,CLM,GAMM] = hinfsyn(P,nmeas,ncont,g);
[K,CL,GAM] = hinfgsyn(P,nmeas,ncont,g);

GAMM
GAM

% bode(KM,K)
% legend('matlab','defsub');

sigma(CLM,CL,ss(g));
legend('matlab','defsub','\gamma');