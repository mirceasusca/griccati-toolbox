%
%   Define nominal plant
%
clearvars
clc;

C = ureal('C',6.e-4,'percentage',[-20,20]); 
C1 = C; C0 = C; C2 = C/10;
Rld = ureal('Rld',15,'percentage',[-20,20]);
L = ureal('L',4.e-5,'percentage',[-20,20]); L0 = L;
E = 12;
Vin = 12;
r_L = 0.01;
rL = r_L;
%
r_C = .2; rC = r_C;
% % 
rC2 = r_C;
rC1 = r_C;
% %
R = .1;
% r_DS = .1;
r_DS = .01; rDS = r_DS;
% % 
r_DS2 = r_DS; rDS2 = r_DS2;
r_DS1 = r_DS; rDS1 = r_DS1;
% 
V_F = .2; VF = V_F;
VF1 = V_F;
VF2 = V_F;
%
% U_ref = 37.5;
U_ref = 24;
% Id = 3.3131;
% mu_d = .623;
mu_d = .5;
% mu_d = .5;
%
Tsc = 1.e-4;
% Ts = 5.e-6;
fPWM=50e3;
TPWM = 1/fPWM;
Ts = 1/fPWM/1000;
% Tsc = 1.e-4;
%
A_on = [-(r_L+r_DS)/L, 0; 0, -1/(r_C+Rld)/C];
B_on = [1/L; 0];
% C_on = eye(2);
% D_on = zeros(2,1);
C_on = [0, 1];
D_on = [0];
%
A_off = [-(r_L+r_DS+r_C*Rld/(r_C+Rld))/L, -Rld/(r_C+Rld)/L; Rld/(r_C+Rld)/C, -1/(r_C+Rld)/C];
B_off = [1/L; 0];
% C_off = eye(2);
% D_off = zeros(2,1);
C_off = [0, 1];
D_off = [0];
%
A_avg = A_on*mu_d + (1-mu_d)*A_off;
B_avg = B_on*mu_d + (1-mu_d)*B_off;
C_avg = C_on*mu_d + (1-mu_d)*C_off;
D_avg = D_on*mu_d + (1-mu_d)*D_off;
%
ny = 1; nu = 1;
%
%   Performance weight
%
s = tf('s');
%
% W1 = (1.414 s + 4)/(2 s + 0.04);
% [Aw1,Bw1,Cw1,Dw1] = tf2ss([14.14, 1600],[40, 16]);
[Aw1,Bw1,Cw1,Dw1] = tf2ss([16.67, 400],[20, 4]);
W1 = ss(Aw1,Bw1,Cw1,Dw1);
%
%   Robustness weight 
%
% W2 = (40*s + 7.2)/(s + 7200);
[Aw2,Bw2,Cw2,Dw2] = tf2ss([40, 2],[1, 2000]);
W2 = ss(Aw2,Bw2,Cw2,Dw2);
%
%   Define augmented plant 
%
%   L = Ln +/-20*Ln
%
delta_L = 0.2;
A = blkdiag(A_avg,Aw1,Aw2); A(3,2) = -1; A(4,2) = 1;
B = zeros(4,4); B(1,3) = -delta_L; B(3,1) = 1; B(3,2) = -1; B(1,4) = B_avg(1,1);
C = [ 0, -Dw1, Cw1,   0;
      0,  Dw2,   0, Cw2;
      A_avg(1,:),0,   0;
      0,   -1,   0,   0];
D = zeros(3,3); D(1,1) = Dw1; D(1,2)=-Dw1; 
D(3,3) = -.1; D(3,4) = B_avg(1,1); 
D(4,1) = 1; D(4,2) = -1;
P = ss(A,B,C,D);

% G = ss(A_avg,B_avg,C_avg,D_avg);
% P = augw(G,W1,[],W2);
% 
%   Optimal controller with hinfsyn
%
opts = hinfsynOptions('Regularize','off')
P = P.NominalValue;
[KM,clpm,gamHinfsynm] = hinfsyn(P,ny,nu);
% [K,clp,gamHinfsyn] = hinfgsyn(P,ny,nu,2.7);
[K,clp,gamHinfsyn] = hinfgsyn(P,ny,nu,2.6988);
%
% [KM,clpm,gamHinfsynm] = hinfsyn(P,ny,nu,[2.5,2.7],opts); % with L
% [K,clp,gamHinfsyn] = hinfgsyn(P,ny,nu,2.7);
gamHinfsynm
gamHinfsyn
subplot(211);
sigma(clpm,clp,ss(gamHinfsynm),ss(gamHinfsyn))
legend('clpm','clp','\gamma_m','\gamma')
% 
% %%
% G = ss(A_avg,1/(1-mu_d)*B_avg,[0,1],0);
% G.InputName = {'mu'};
% G.OutputName = {'uC'};
% K.InputName = {'uC'};
% K.OutputName = {'mu'};
% KM.InputName = {'uC'};
% KM.OutputName = {'mu'};
% CL = feedback(G,K);
% CLM = feedback(G,KM);
% 
% subplot(212);
% step(CLM,CL)
% legend('CLM','CL')
