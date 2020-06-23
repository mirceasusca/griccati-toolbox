A = [-4.019, 5.12, 0, 0, -2.082, 0, 0, 0, 0.87;
    -0.346, 0.986, 0, 0, -2.34, 0, 0, 0, 0.97;
    -7.909, 15.407, -4.069, 0, -6.45, 0, 0, 0, 2.68;
    -21.816, 35.606, -0.339, -3.87, -17.8, 0, 0, 0, 7.39;
    -60.196, 98.188, -7.907, 0.34, -53.008, 0, 0, 0, 20.4;
    0, 0, 0, 0, 94.0, -147.2, 0, 53.2, 0;
    0, 0, 0, 0, 0, 94.0, -147.2, 0, 0;
    0, 0, 0, 0, 0, 12.8, 0, -31.6, 0;
    0, 0, 0, 0, 12.8, 0, 0, 18.8, -31.6];

B = [0.01, 0.003, 0.009, 0.024, 0.068, 0, 0, 0, 0;
    -0.011, -0.021, -0.059, -0.162, -0.445, 0, 0, 0, 0;
    -0.151, 0, 0, 0, 0, 0, 0, 0, 0]';

Q = eye(size(A)); R = eye(size(B,2)); L = zeros(size(B));

sigma = create_popov_triplet(A,B,Q,L,R);

[F,Xk,info] = gricsv(sigma);

Xk = gricir(sigma,Xk);

% X0-X_ideal
% X1-X_ideal

%% Example 2 [17, Example 3], [16, Example 6.15]
% this example illustrates a linear-quadratic control problem
A = [0.9512, 0; 0. 0.9048];
B = [4.877, 4.877; -1.1895, 3.569];
R = [1/3, 0; 0, 3];
Q = [0.005, 0; 0, 0.02];
L = zeros(size(B));

sigma = create_popov_triplet(A,B,Q,L,R,true);
[F,Xk,info] = gricsv(sigma);

Xk = gricir(sigma,Xk);
