function [F,X,info] = lqe_df(A,Bw,C,Dw,Q,R,N,bal)
%LQE_DF Solves the linear-quadratic estimator (LQE) problem/Kalman filter 
% for the continuous-time state-space system:
%   dx/dt = Ax + Bu + Bw*w
%       y = Cx + Du + Dw*w + v
% with known inputs u, white process noise w, and white measurement noise v
% satisfying:
%  E(w) = E(v) = 0, E(ww') = Q, E(vv') = R, E(wv') = N.
% It constructs a state estimate xhat(t) which minimizes the steady-state
% error covariance:
%  P = lim E({x-xhat}{x-xhat}').
%
% info.e = eig(A+F*C);
%
% The estimator equation is:
%  dxhat/dt = A*xhat + B*u - F(y-C*xhat-D*u).
%
% [F,X,INFO] = LQE_DF(A,Bw,C,Dw,Q,R)
% [F,X,INFO] = LQE_DF(A,Bw,C,Dw,Q,R,N,'balance')
% [F,X,INFO] = LQE_DF(A,Bw,C,Dw,Q,R,N,'nobalance')
%

if nargin == 6
    N = zeros(size(Q*Dw'));
    bal = 'balance';
elseif nargin == 7
    bal = 'balance';
end

% m2 = size(G,2);
% p2 = size(H,2);

Q_bar = Bw*Q*Bw';
R_bar = R+Dw*N+N'*Dw'+Dw*Q*Dw';
N_bar = Bw*(Q*Dw'+N);
sigma1 = create_popov_triplet(A',C',Q_bar,N_bar,R_bar);

[F,X,info]=gricsv(sigma1,bal);
F = F';
info.e = eig(A+F*C);

end
