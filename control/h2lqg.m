function [K,CL,GAM,INFO] = h2lqg(A,B,C,D,Qx,Ru,Bw,Qw,Rv,discr,bal)
%H2LQG Solves the LQG problem using the H2 framework.
%
%      [      A  | Bw*sqrt(Qw)  0  |       B ]
%      [ ----------------------------------- ]
%  P = [ sqrt(Qx)|      0       0  |       0 ]
%      [      0  |      0       0  | sqrt(Ru)]
%      [ ----------------------------------- ]
%      [      C  |      0  sqrt(Rv)|       D ]
%
%  The plant model is described by:
% dx = Ax + Bu + Bw*w;
%  y = Cx + Du + v,
% where dx means dx/dt or x[k+1] if the plant is continuous or
% discrete-time, respectively.
%
%  Noise signal assumptions:
% E(w) = E(v) = 0, E(ww') = Qw, E(vv') = Rv, E(wv') = 0.
% 
%  [K,CL,GAM,INFO] = h2lqg(A,B,C,D,Qx,Ru,Bw,Qw,Rv)
%  [K,CL,GAM,INFO] = h2lqg(A,B,C,D,Qx,Ru,Bw,Qw,Rv,discr)
%  [K,CL,GAM,INFO] = h2lqg(A,B,C,D,Qx,Ru,Bw,Qw,Rv,discr,bal)
% 

if nargin == 9
    discr = false;
    bal = 'balance';
elseif nargin == 10
    bal = 'balance';
end

[P,p2,m2] = create_lqg_h2_gen_plant(A,B,C,D,Qx,Ru,Bw,Qw,Rv,discr);
[K,CL,GAM,INFO] = h2gsyn(P,p2,m2,bal);
end

function [P,p2,m2] = create_lqg_h2_gen_plant(A,B,C,D,Qx,Ru,Bw,Qw,Rv,discr)
%CREATE_LQG_H2_GEN_PLANT Creates LQG/H2 generalized plant P.
%

n = size(A,1);
%
nw = size(Qw,1);
nv = size(Rv,1);
%
nx = size(Qx,1);
ny = size(C,1);
nu = size(B,2);
%
p2 = ny;
m2 = nu;

B1 = [Bw*real(sqrtm(Qw)),zeros(n,nv)];
C1 = [real(sqrtm(Qx)); zeros(nu,nx)];
D11 = zeros(nx+nu,nw+nv);
D12 = [zeros(nx,nu);sqrtm(Ru)];
D21 = [zeros(ny,nw),sqrtm(Rv)];
D22 = D;

P = ss(A,[B1,B],[C1;C],[D11,D12;D21,D22],discr*1);
end