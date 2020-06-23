function [F,X,info] = dlqr_df(A,B,Q,R,L,bal)
%DLQR_DF Solves the linear-quadratic regulator (LQR) problem for 
% discrete-time systems. Computes the state-feedback law u = F*x which 
% minimizes the quadratic cost function:
%  J = sum_{n=0}^{inf} (x'Qx + u'Ru + 2x'Lu).
%
% info.e = eig(A+B*F)
%
% [F,X,INFO] = DLQR_DF(A,B,Q,R)
% [F,X,INFO] = DLQR_DF(A,B,Q,R,L,'balance')
%

if nargin == 4
    L = zeros(size(Q,1),size(R,2));
    bal = 'balance';
elseif nargin == 5
    bal = 'balance';
end

sigma1 = create_popov_triplet(A,B,Q,L,R,true);

% check solvability conditions

[F,X,info]=gricsv(sigma1,bal);
info.e = eig(A+B*F);

end
