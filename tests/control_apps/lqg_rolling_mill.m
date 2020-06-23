%%
clearvars
clc

%%
% Hydraulic actuator (with input "u-x")
Hx = tf(2.4e8,[1 72 90^2],'inputname','u-x');

% Input thickness/hardness disturbance model
Fix = tf(1e4,[1 0.05],'inputn','w-ix');

% Rolling eccentricity model
Fex = tf([3e4 0],[1 0.125 6^2],'inputn','w-ex');

% Gain from force to thickness gap
gx = 1e-6;

% I/O map from inputs to forces f1 and f2
Px = append([ss(Hx) Fex],Fix);

% Add static gain from f1,f2 to outputs "x-gap" and "x-force" 
Px = [-gx gx;1 1] * Px;

% Give names to the outputs:
set(Px,'outputn',{'x-gap' 'x-force'});

lpf = tf(30,[1 30]);

% Connect low-pass filter to first output of Px
Pxdes = append(lpf,1) * Px;
set(Pxdes,'outputn',{'x-gap*' 'x-force'});

% Design the state-feedback gain using LQRY and q=1, r=1e-4
kx = lqry(Pxdes(1,1),1,1e-4);

[estx,L] = kalman(Pxdes(2,:),eye(2),1000);

Regx = lqgreg(estx,kx);

% h = bodeplot(Regx,{0.1 1000})
% setoptions(h,'PhaseMatching','on')

clx = feedback(Px,Regx,1,2,+1);       % Note: +1 for positive feedback

% h = bodeplot(Px(1,2:3),'--',clx(1,2:3),'-',{0.1 100})
% setoptions(h,'PhaseMatching','on')

dt = 0.01;
t = 0:dt:50;   % time samples

% Generate unit-covariance driving noise wx = [w-ex;w-ix].
% Equivalent discrete covariance is 1/dt
wx = sqrt(1/dt) * randn(2,length(t));

lsim(Px(1,2:3),':',clx(1,2:3),'-',wx,t);

%%
% Specify model components
Hy = tf(7.8e8,[1 71 88^2],'inputn','u-y');
Fiy = tf(2e4,[1 0.05],'inputn','w-iy');
Fey = tf([1e5 0],[1 0.19 9.4^2],'inputn','w-ey');
gy = 0.5e-6; % force-to-gap gain

% Build open-loop model
Py = append([ss(Hy) Fey],Fiy);
Py = [-gy gy;1 1] * Py;
set(Py,'outputn',{'y-gap' 'y-force'});

% State-feedback gain design
Pydes = append(lpf,1) * Py; % Add low-freq. weigthing
set(Pydes,'outputn',{'y-gap*' 'y-force'});
ky = lqry(Pydes(1,1),1,1e-4);

% Kalman estimator design
esty = kalman(Pydes(2,:),eye(2),1e3);

% Form SISO LQG regulator for y-axis and close the loop
Regy = lqgreg(esty,ky);
cly = feedback(Py,Regy,1,2,+1);

dt = 0.01;
t = 0:dt:50;
wy = sqrt(1/dt) * randn(2,length(t));

lsim(Py(1,2:3),':',cly(1,2:3),'-',wy,t);

%% Model couping
P = append(Px,Py);
P = P([1 3 2 4],[1 4 2 3 5 6]);

gxy = 0.1; gyx = 0.4;
CCmat = [eye(2) [0 gyx*gx;gxy*gy 0] ; zeros(2) [1 -gyx;-gxy 1]];
Pc = CCmat * P;
Pc.outputname = P.outputname;

%%
feedin = 1:2 % first two inputs of Pc are the commands
feedout = 3:4 % last two outputs of Pc are the measurements
cl = feedback(Pc,append(Regx,Regy),feedin,feedout,+1)

%%
Pdes = append(lpf,lpf,eye(2)) * Pc;
Pdes.outputn = Pc.outputn;

Qx = eye(2);
Ru = 1e-4*eye(2);
k = lqry(Pdes(1:2,1:2),Qx,Ru);      % LQ gain

% own lqr function
Pdes_lqr = Pdes(1:2,1:2);
[Qxx,Lxxu,Rxu]=lqry_to_lqr_mat(Qx,zeros(size(Qx,2),size(Ru,1)),Ru,Pdes_lqr.C,Pdes_lqr.D);
K2 = lqr_df(Pdes_lqr.A,Pdes_lqr.B,Qxx,Rxu,Lxxu);

Qw = eye(4);
Rv = 1e3*eye(2);
[est,L] = kalman(Pdes(3:4,:),Qw,Rv);     % Kalman estimator

% own lqe function
P_lqe = Pdes(3:4,:);
L2 = lqe_df(P_lqe.a,P_lqe.b(:,3:end),P_lqe.c,P_lqe.d(:,3:end),Qw,Rv);

A = Pdes.a;
B = Pdes.b(:,1:2);
C = Pdes.c(3:4,:);
D = Pdes.d(3:4,1:2);
K_lqg = ss(A+B*K2+L2*C,-L2,K2,0);

%%
RegMIMO = lqgreg(est,k);      % form MIMO LQG regulator

% sigma(K_lqg,RegMIMO)

% Form the closed-loop model
cl1 = feedback(Pc,RegMIMO,1:2,3:4,+1);
cl2 = feedback(Pc,K_lqg,1:2,3:4,+1);

wxy = [wx ; wy];

figure
% Simulate with lsim using same noise inputs
subplot(121)
lsim(Pc(1:2,3:6),':',cl1(1:2,3:6),'-',wxy,t);
subplot(122)
lsim(Pc(1:2,3:6),':',cl1(1:2,3:6),'-',wxy,t);