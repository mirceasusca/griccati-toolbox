%% Parameters
xi = 0.05;
alpha = [0.09877, -0.309, -0.891, 0.5878, 0.7071, -0.8091];
w = [1, 4, 9, 16, 25, 36];

%% Beam model
G = tf(alpha(1)^2*[1,0],[1, 2*xi*w(1), w(1)^2]) + ...
    tf(alpha(2)^2*[1,0],[1, 2*xi*w(2), w(2)^2]) + ...
    tf(alpha(3)^2*[1,0],[1, 2*xi*w(3), w(3)^2]) + ...
    tf(alpha(4)^2*[1,0],[1, 2*xi*w(4), w(4)^2]) + ...
    tf(alpha(5)^2*[1,0],[1, 2*xi*w(5), w(5)^2]) + ...
    tf(alpha(6)^2*[1,0],[1, 2*xi*w(6), w(6)^2]);

G.InputName = 'uG';  G.OutputName = 'y';

[A,B,C,D] = ssdata(G);

%% LQR
Qx = C'*C;
Ru = 1e-3;

[Kr,Xr,info_r] = lqr_df(A,B,Qx,Ru);

%% LQE
Qw = 1;
Rv = 1e-2;

% dx/dt = Ax + B(u+w); y = Cx + D(u+w) + v
[Kf,Xf,info_f] = lqe_df(A,B,C,D,Qw,Rv);

%% LQG controller using GRICSV
K = ss(A+B*Kr+Kf*C,-Kf,Kr,0);

%% LQG controller using H2 GRICSV
[KLQG,CL,GAM,INFO] = h2lqg(A,B,C,D,Qx,Ru,B,Qw,Rv);

%% LQG controller MATLAB
M = [C D;zeros(1,12) 1];  % [y;u] = M * [x;u]
QWV = blkdiag(B*B',1e-2);
QXU = M'*diag([1 1e-3])*M;
G = ss(G);
[CLQG,info_lqg] = lqg(ss(G),QXU,QWV);

%%
step(K,KLQG,CLQG)

%%
S1 = sumblk('yn = y + n');
S2 = sumblk('uG = u + d + r');
%
K.InputName = {'yn'};
K.OutputName = {'u'};
%
CL0 = connect(G,K,S1,S2,{'d','n','r'},{'y','u'});
G0 = CL0({'y'},{'r'});

%%
T = 1e-4;
t = 0:T:60;
w = 8.98;
u = sin(w*t);
lsim(G,G0,u,t)
legend('G','CL')

%%
% subplot(211)
% impulse(G,G0)
%
% subplot(212)
% bodemag(CLQG,K)

