%
% Q = C'*Qt*C;
% G = B/R*B'

%% Example 1 [27, Example 1]
A = [0,1;0,0]; B = [0;1]; R = 1; Q = [1,0;0,2]; L = [0;0];
sigma = create_popov_triplet(A,B,Q,L,R);

[M,N] = create_hamiltonian_pencil(sigma);
tol = size(N,1)*size(N,2)*eps(max(1,norm(N,1)));
[F,X,info] = gricsv(sigma);

C = abs(X - [2,1;1,2]) < tol';
assert(all(C(:)) == 1)

X1 = care(A,B,Q,R,L);
F1 = -R\B'*X;

% the solution from gricsv is identical to the solution from care
eig(A+B*F) % ideally {-1,-1}; 
eig(A+B*F1) % ideally {-1,-1}; 
% assert()

%% Example 2 [27, Example 2]
A = [4,3;-9/2,-7/2]; B = [1;-1]; R = 1; Q = [9,6;6,4]; L = zeros(size(B));
sigma = create_popov_triplet(A,B,Q,L,R);

[M,N] = create_hamiltonian_pencil(sigma);
tol = size(N,1)*size(N,2)*eps(max(1,norm(N,1)));
[F,X,info] = gricsv(sigma);

X_ideal = [
    9*(1+sqrt(2)), 6*(1+sqrt(2));
    6*(1+sqrt(2)), 4*(1+sqrt(2))];

C = X-X_ideal;

assert(all((abs(C(:))) < tol))

eig(A+B*F) % should be {-sqrt(2), -1/2}
X
X_ideal

sigma_new = create_popov_triplet(A+B*F,B,Q,L,R);
[F2,X2,info2] = gricsv(sigma_new)
F2
X2

%% Example 3 [7]
A = [
    0, 1, 0, 0;
    0, -1.89, 0.39, -5.53;
    0, -0.034, -2.98, 2.43;
    0.034, -0.0011, -0.99, -0.21];

B = [0, 0;
    0.36, -1.6;
    -0.95, -0.032;
    0.03, 0];

Q = [2.313, 2.727, 0.688, 0.023;
    2.727, 4.271, 1.148, 0.323;
    0.688, 1.148, 0.313, 0.102;
    0.023, 0.323, 0.102, 0.083];

R = eye(2);

L = zeros(size(B));

sigma = create_popov_triplet(A,B,Q,L,R);
    
eig(Q) % in this example, Q has a small eigenvalue of order O(1e-3), which
% may reflect a small perturbation in the input data. The computed 
% stabilizing solution is nevertheless positive definite with eigenvalues
% greater than one (??).

[F,X,info] = gricsv(sigma);
% info
% eig(X)

% X1 = care(A,B,Q,R,L);
% eig(X1)

eig(A+B*F)

%% Example 4 
% Mathematical model of a binary distillation column with condenser,
% reboiler and nine plates
d0 = [-0.991,-1.051,-1.118,-1.548,-1.640,-1.721,-1.823,-1.943];
d1 = [0.529,0.596,0.596,0.718,0.799,0.901,1.021];
dm1 = [0.522,0.522,0.522,0.922,0.922,0.922,0.922];
A = diag(d0) + diag(d1,1) + diag(dm1,-1);

B = 1e-3*[3.84,4,37.6,3.08,2.36,2.88,3.08,3;
    -2.88,-3.04,-2.80,-2.32,-3.32,-3.82,-4.12,-3.96]';

R = eye(2);
Q = [
    1,0,0,0,0.5,0,0,0.1;
    0,1,0,0,0.1,0,0,0;
    0,0,1,0,0,0.5,9,0;
    0,0,0,1,0,0,0,0;
    0.5,0.1,0,0,0.1,0,0,0;
    0,0,0.5,0,0,0.1,0,0;
    0,0,9,0,0,0,0.1,0; % 0,0,0,0,0,0,0.1,0
    0.1,0,0,0,0,0,0,0.1];

L = eye(size(B));

% Q = (Q+Q')/2; % symmetrize matrix
sigma = create_popov_triplet(A,B,Q,L,R);
[F,X,info] = gricsv(sigma);
info

% note that Q is indefinite and the computed stabilizing solution is also
% indefinite
eig(Q)
eig(X)
eig(A+B*F)

[s1,s2,s3] = validate_gars(sigma,F,X,info.V,info.r);
norm(s1)
norm(s2)
norm(s3)
% % 
X1 = care(A,B,Q,R,L);
F1 = -R\(B'*X+L');
[s1,s2,s3] = validate_gars(sigma,F1,X1);
norm(s1)
norm(s2)
norm(s3)

%% Example 5 [34]
% This is the data for a 9th-order continuous state space model of a 
% tubular ammonia reactor. It should be noted that the underlying model 
% includes a disturbance term which is neglected in this context.
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

[M,N] = create_hamiltonian_pencil(sigma);
cond(N)
norm(N)
norm(X)
[F,X,info] = gricsv(sigma);
info

eig(X)
eig(A)
eig(A+B*F)

%% Example 7 -- parameter-dependent problems of fixed size
e = 1e-6;
A = [1,0;0,-2]; B = [e;0];
R = 1; C = [1,1]; Qt = 1;
Q = C'*Qt*C;
L = zeros(size(B));

Sigma = create_popov_triplet(A,B,Q,L,R);
[F,X,info] = gricsv(Sigma);

X_sol = @(e)[(1+sqrt(1+e^2))/e^2,1/(2+sqrt(1+e^2));
    1/(2+sqrt(1+e^2)),1/4*(1-e^2/(2+sqrt(1+e^2))^2)];
X-X_sol(e)

%% Example 8
e = 1;
% e = 1e-8;
A = [-0.1,0;0,-0.02]; B = [0.1,0;0.001,0.01];
R = [1+e,1;1,1]; C = [10,100]; Qt = 1;
Q = C'*Qt*C;
L = zeros(size(B));

Sigma = create_popov_triplet(A,B,Q,L,R);
[F,X,info] = gricsv(Sigma);
cond(X)

%% Example 9
e = 1e-3;
A = [0,e;0,0];
B = [0;1];
R = 1;
Q = eye(2);
L = zeros(size(B));

Sigma = create_popov_triplet(A,B,Q,L,R);
[F,X,info] = gricsv(Sigma);

X_sol = @(e) [sqrt(1+2*e)/e,1;1,sqrt(1+2*e)];
eig(X_sol(e))
eig(X)
eig(X_sol(e))-eig(X)

%% Example 10
e = 1;
A = [e+1,1;1,e+1];
B = eye(2);
R = eye(2); % G = eye(2); G = B/R*B'
Q = diag([e^2,e^2]);
L = zeros(size(B));

x11 = 1/2*(2*(e+1)+sqrt(2*(e+1)^2+2) + sqrt(2)*e);
x22 = x11;
x12 = x11/(x11-(e+1));
x21 = x12;
X_sol = [x11,x12;x21,x22];

Sigma = create_popov_triplet(A,B,Q,L,R);
[F,X,info] = gricsv(Sigma);

%% Example 11 -- example arising in Hinf control problem
e = 2;
A = [3-e,1;4,2-e];
B = [1;1];
R = 1;
Q = [4*e-11,2*e-5;2*e-5,2*e-2];
L = zeros(size(B));

Sigma = create_popov_triplet(A,B,Q,L,R);
[F,X,info] = gricsv(Sigma);
X % -- invariant

%% Example 12
e = 2;

v = [1,1,1]';
V = eye(3) - 2/3*(v*v');
A0 = e*diag([1,2,3]);
Q0 = diag([1/e,1,e]);
A = V*A0*V;
%
B = eye(3);
R = e*eye(3);
%
Q = V*Q0*V;
L = zeros(size(B));

Sigma = create_popov_triplet(A,B,Q,L,R);
[F,X,info] = gricsv(Sigma);

X_sol = V*diag([e^2+sqrt(e^4+1);
    2*e^2+sqrt(4*e^4+e);
    3*e^2+sqrt(9*e^4+e^2)])*V;

X-X_sol

%% Example 13 -- magnetic tape control problem
e = 0.1;
A = [0,0.4,0,0;
    0,0,0.345,0;
    0,-0.524/e,-0.465/e,0.262/e;
    0,0,0,-1/e];
B = [0;0;0;1/e];
Q = diag([1,0,1,0]);
R = 1;
L = zeros(size(B));

Sigma = create_popov_triplet(A,B,Q,L,R);
[F,X,info] = gricsv(Sigma);
X
cond(X)

%% Example 14
e = 1e-6;
A = [-e,1,0,0;
    -1,-e,0,0;
    0,0,e,1;
    0,0,-1,e];
B = [1;1;1;1];
R = 1;
C = [1,1,1,1];
Qt = 1;
Q = C'*Qt*C;
L = zeros(size(B));

Sigma = create_popov_triplet(A,B,Q,L,R);
[F,X,info] = gricsv(Sigma);
X
cond(X)

%% Example 15 -- examples of scalable size without parameters
