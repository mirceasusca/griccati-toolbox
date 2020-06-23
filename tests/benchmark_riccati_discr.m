%% Example -- On computing the stabilizing solution
% of the discrete-time riccati equation - V. Ionescu, M. Weiss
A = [0,1;0,-1];
B = [1,0;2,1];
Q = [-4,-4;-4,7]/11;
L = [3,1;-1,7];
R = [9,3;3,1];
Sigma = create_popov_triplet(A,B,Q,L,R,true);

[F,X,info] = gricsv(Sigma);

RCD = @(X) A'*X*A-X+Q-(L+A'*X*B)/(R+B'*X*B)*(L'+B'*X*A);
norm(RCD(X))

%% Example 1 [17,Example 2]
% this is an example of stabilizable-detectable, but
% uncontrollable-unobservable data
A = [4,3;-9/2,-7/2];
B = [1;-1];
R = 1;
Q = [9,6;6,4];
L = zeros(size(B));

sigma = create_popov_triplet(A,B,Q,L,R,true);

[F1,X1,info] = gricsv(sigma)
X1
X_analytic = (1+sqrt(5))/2*[9,6;6,4]
[S1,S2,S3] = validate_gars(sigma,F1,X1)

[M,N] = create_hamiltonian_pencil(sigma);
X2 = gdare(N,M,size(A,1))

%% Example 2 [17, Example 3], [16, Example 6.15]
% this example illustrates a linear-quadratic control problem
A = [0.9512, 0; 0. 0.9048];
B = [4.877, 4.877; -1.1895, 3.569];
R = [1/3, 0; 0, 3];
Q = [0.005, 0; 0, 0.02];
L = zeros(size(B));

sigma = create_popov_triplet(A,B,Q,L,R,true);

[F1,X1,info] = gricsv(sigma)
X1
[M,N] = create_hamiltonian_pencil(sigma);
X2 = gdare(N,M,size(A,1))