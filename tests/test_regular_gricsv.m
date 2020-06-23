%%
rng(8);
n = 30;
m = 5;
% 
A = rand(n)+diag(1*rand(1,n-2),2);
% A = rand(n);
B = rand(n,m);

Q = rand(n);
Q = (Q+Q')/2 + 10*eye(n);
% L = zeros(n,m);
L = 0.01*rand(n,m);
R = eye(m);

discr = true;
margin = 0.95;
sigma = create_popov_triplet(A,B,Q,L,R,discr,margin);

[Fgricsv,Flqr,Fgare]=compare_gricsv_with_lqr(sigma);

%%
rng(8);
n = 300; % (38,2), (280,10), (280,11), (400,35), (40,3) (36,1)
m = 60;
%
% A = rand(n)+diag(100*rand(1,n-2),2);
A = rand(n);
B = rand(n,m);

Q = rand(n);
Q = (Q+Q')/2 + 10*eye(n);
% L = zeros(n,m);
L = 0.01*rand(n,m);
R = eye(m);

discr = false;
margin = -0.2;
sigma = create_popov_triplet(A,B,Q,L,R,discr,margin);

[Fgricsv,Flqr,Fgare] = compare_gricsv_with_lqr(sigma);
