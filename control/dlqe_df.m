function [F,X,info] = dlqe_df(A,G,C,H,Q,R,N,bal)
%DLQE_DF Solves the linear-quadratic estimator (LQE) problem/Kalman filter 
% for the discrete-time state-space system:
%   x[n+1] = Ax[n] + Bu[n] + Gw[n]
%     y[n] = Cx[n] + Du[n] + Hw[n] + v[n]
% with known inputs u, white process noise w, and white measurement noise v
% satisfying:
%  E(w) = E(v) = 0, E(ww') = Q, E(vv') = R, E(wv') = N.
% It constructs a state estimate xhat[k] which minimizes the steady-state
% error covariance:
%  P = lim E({x-xhat}{x-xhat}').
%
% info.e = eig(A+F*C);
%
% The estimator equation is:
%  xhat[n+1] = A*xhat[n] + B*u[n] - F(y[n]-C*xhat[n]-D*u[n]).
%
% [F,X,INFO] = DLQE_DF(A,G,C,H,Q,R)
% [F,X,INFO] = DLQE_DF(A,G,C,H,Q,R,N,'balance')
% [F,X,INFO] = DLQE_DF(A,G,C,H,Q,R,N,'nobalance')
%

if nargin == 6
    N = zeros(size(Q*H'));
    bal = 'balance';
elseif nargin == 7
    bal = 'balance';
end

% m2 = size(G,2);
% p2 = size(H,2);

Q_bar = G*Q*G';
R_bar = R+H*N+N'*H'+H*Q*H';
N_bar = G*(Q*H'+N);
sigma1 = create_popov_triplet(A',C',Q_bar,N_bar,R_bar,true);

% check solvability conditions

[F,X,info]=gricsv(sigma1,bal);
F = F';
info.e = eig(A+F*C);

end
