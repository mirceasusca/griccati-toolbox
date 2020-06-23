function [Q1,L1,R1] = lqry_to_lqr_mat(Q,L,R,C,D)
%LQRY_TO_LQR_MAT
%
% int_0^{\infty} y'*Q*y + u'*R*u+2*y'*L*u
% int_0^{\infty} x'*Q1*x+u'*R1'u+2*x'*L1'*u
%
% y = Cx+Du.
%

Q1 = C'*Q*C;
L1 = C'*Q*D+C'*L;
R1 = D'*Q*D+L'*D+D'*L+R;
end