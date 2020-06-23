A = [-0, 1; 0, -2];
B = [1e-7;0];

Q = eye(2);
L = zeros(size(B));
R = 1;

F = lqr_df(A,B,Q,R,L)
eig(A)
eig(A+B*F)