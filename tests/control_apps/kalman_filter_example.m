%% 1) create ss system
A = [0 1 0;0 0 1;1 0 0];    
B = [0.3 1;0 1;-0.3 0.9];
C = [1.9 1.3 1];  
D = [0.53 -0.61];
sys = ss(A,B,C,D);

%% 2) define state and input weights -- apply lqr
nx = 3;
Qx = 0.1*eye(nx);
Ru = [1, 0; 0, 2];

[KR1,S1,E1] = lqr(A,B,Qx,Ru)
[KR2,S2,info_r] = lqr_df(A,B,Qx,Ru)

%% 3) define noise covariance data and weighting matrices -- apply lqe
G = B(:,2);
H = D(:,2);
Qw = 1; Rv = 1; Nwv = 0;

[~,K1,P1] = kalman(sys,Qw,Rv,Nwv,1,1)
[K2,P2,info_e] = lqe_df(A,G,C,H,Qw,Rv)
