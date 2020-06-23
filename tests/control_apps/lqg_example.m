%% 1) create ss system
A = [0 1 0;0 0 1;1 0 0];    
B = [0.3 1;0 1;-0.3 0.9];
C = [1.9 1.3 1];  
D = [0.53 -0.61];
% D = [0.53 0];
sys = ss(A,B,C,D);

pole(sys)

%% 2) define noise covariance data and weighting matrices
nx = 3; % num. states
ny = 1; % num. outputs
Qn = [4 2 0; 2 1 0; 0 0 1];
Rn = 0.7;
R = [1, 0; 0, 2];
QXU = blkdiag(0.1*eye(nx),R);
QWV = blkdiag(Qn,Rn);
QI = eye(ny);

[KEST,K1,P1] = kalman(sys,1,1,0,1,1)
% [K2,P2,info] = lqe_df(A,C,B(:,2),1,1)
[K2,P2,info] = lqe_df(A,B(:,2),C,D(:,2),1,1)

%%
KLQG = lqg(sys,QXU,QWV);
KLQG1 = lqg(sys,QXU,QWV,QI,'1dof');
KLQG2 = lqg(sys,QXU,QWV,QI,'2dof');

%%
L = series(KLQG2,sys);
CL = feedback(L,1,2,1,+1);
CL = CL(1,1);
pole(CL)
step(CL);