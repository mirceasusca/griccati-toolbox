rng(8); % n = 15, m = 6, balancing has an order of magnitude effect
% n = 45, m = 244; lqr blows
n = 23;
m = 45;
%
A = rand(n)+diag(800*rand(1,n-2),2);
% A = rand(n);
% A = randn(n);
B = rand(n,m);

Q = 1.1*eye(n);
L = zeros(n,m);
R = rand(m);
R = 3*eye(m) + (R+R')/2;

discr = true;
margin = 0.95;
sigma = create_popov_triplet(A,B,Q,L,R,discr,margin);

[Fgricsv,Flqr,Fgare]=compare_gricsv_with_lqr(sigma);
[V,W,J,F,X] = kpys(sigma);
% V

disp('');
[s1,s2,s3] = validate_gars(sigma,F,X);
disp('>> Norms of KPYS experiment');
disp([norm(s1),norm(s2),norm(s3)]);
disp(['Stability of KPYS solution: ',num2str(all(eig(A+B*F) < 1))]);

%%
rng(8);
n = 30; % (38,2), (280,10)w, (280,11)m, (400,35)best, (40,3)b (36,1)b
m = 4;
%
A = rand(n)+diag(100*rand(1,n-2),2);
B = rand(n,m);

Q = eye(n);
L = zeros(n,m);
R = rand(m);
R = eye(m) + (R+R')/2;

discr = false;
margin = -0.2;
sigma = create_popov_triplet(A,B,Q,L,R,discr,margin);

[Fgricsv,Flqr,Fgare]=compare_gricsv_with_lqr(sigma);
[V,W,J,F,X] = kpys(sigma);
% V

disp('');
[s1,s2,s3] = validate_gars(sigma,F,X);
disp('>> Norms of KPYS experiment');
disp([norm(s1),norm(s2),norm(s3)]);
disp(['Stability of KPYS solution: ',num2str(all(eig(A+B*F) < 0))]);