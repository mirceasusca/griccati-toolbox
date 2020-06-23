function [F,X,info] = lqr_df(A,B,Q,R,L,bal)
%LQR_DF Solves the linear-quadratic regulator (LQR) problem for 
% continuous-time systems. Computes the state-feedback law 
% u = F*x which minimizes the quadratic cost function:
%  J = int_0^{inf} (x'Qx + u'Ru + 2x'Lu) dt.
%
% info.e = eig(A+B*F)
%
% [F,X,INFO] = LQR_DF(A,B,Q,R)
% [F,X,INFO] = LQR_DF(A,B,Q,R,L,'balance')
%

if nargin == 4
    L = zeros(size(Q,1),size(R,2));
    bal = 'balance';
elseif nargin == 5
    bal = 'balance';
end

sigma1 = create_popov_triplet(A,B,Q,L,R);

[F,X,info]=gricsv(sigma1,bal);
info.e = eig(A+B*F);

end
