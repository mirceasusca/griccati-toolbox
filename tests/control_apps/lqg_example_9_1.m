G = nd2sys([-2,1],[50,15,1],3);
int = nd2sys(1,[1,0.0001]);
%
Gs = mmult(G,int);
[A,B,C,D] = unpck(Gs);
Gs = ss(A,B,C,D);
Q = 1.0*(C'*C);
R = 1;
%
Kx = lqr(A,B,Q,R);
Kx2 = lqr_df(A,B,Q,R);
Kx+Kx2
%
Bnoise = eye(size(A));
W = eye(size(A));
V = 1;
%
Ke = lqe(A,Bnoise,C,W,V);
Ke2 = lqe_df(A,Bnoise,C,zeros(size(C,1),size(Bnoise,2)),W,V);
Ke+Ke2
%
[Ac,Bc,Cc,Dc] = reg(A,B,C,D,Kx,Ke);
Klqg1 = ss(Ac,Bc,-Cc,Dc);
Ks = pck(Ac,Bc,Cc,Dc);
Klqg = mmult(Ks,int);
%
[Ak,Bk,Ck,Dk] = unpck(Klqg);
K = ss(Ak,Bk,Ck,Dk);

[Klqg2,CL,GAM,INFO] = h2lqg(A,B,C,D,Q,R,Bnoise,W,V);

subplot(211)
bode(Klqg1,Klqg2)

subplot(212)
step(Klqg1,Klqg2)