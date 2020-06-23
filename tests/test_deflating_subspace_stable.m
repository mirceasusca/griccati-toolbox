% rng(4);
% Q9 = qr(rand(9));
% Q10 = qr(rand(10));
M = diag([0,0,-1,2,3]);
%
M = [zeros(5,1),M]; % right kron. str L(0,1)
%
M = [zeros(1,size(M,2));M]; % right kron str L(1,2)
M = [zeros(6,2),M];
M(1,1) = 1;
%
M = [M,zeros(6,2)]; % left kron str LT(1,2)
M = [M;zeros(3,10)]; M(7,9) = 1;M(8,10) = 1;
%
M = [M;zeros(1,size(M+1,2))]; % left kron str LT(0,1)
%
% M = Q9*M*Q10;
%
N = diag([5,4,3,2,1]); % right kron. str L(0,1)
%
N = [zeros(5,1),N];
%
N = [zeros(1,size(N,2));N]; % right kron str L(1,2)
N = [zeros(6,2),N];
N(1,2) = 1;
%
N = [N,zeros(6,2)]; % left kron str LT(1,2)
N = [N;zeros(3,10)]; N(8,9) = 1;N(9,10) = 1;
%
N = [N;zeros(1,size(N+1,2))]; % left kron str LT(0,1)

% N = Q9*N*Q10;
parse_gklf(N,M)

[V,info]=mpdefsub(M,N);
S = info.S;
M*V*S-N*V
V
eig(S)
[r,nu_r,nu_l] = nrankp(N,M)

%%
n = 14;
rng(5);
M = rand(n);
N = rand(n);
eig(N,M)
parse_gklf(N,M)
[V,info]=mpdefsub(M,N);
V
S = info.S;
M*V*S-N*V
eig(S)
[r,nu_r,nu_l] = nrankp(N,M)

%%

% num = 10;
M = eye(4);
M(4,4) = 0;
% N = [-2,0,0;0,4,1;0,-1,4];
N = [-2,0,0,0;0,4,1,0;0,0,4,0;0,0,0,5];

parse_gklf(N,M)
[V,info]=mpdefsub(M,N);
S = info.S;
M*V*S-N*V
V
eig(S)
parse_gklf(N,M)
[r,nu_r,nu_l] = nrankp(N,M)

%%
M = [1,0, 0, 0, 0,0,0,0, 0;
     0,0, 0, 0, 1,0,0,0, 0;
     0,0, 0, 0, 0,1,0,0, 0;
     0,0, 0, 0, 0,0,1,0, 0;
     0,0, 0, 0, 0,0,0,0, 3];
 N =[0,1, 0, 0, 0,0,0,0, 0;
     0,0, 0, 0, 0,1,0,0, 0;
     0,0, 0, 0, 0,0,1,0, 0;
     0,0, 0, 0, 0,0,0,1, 0;
     0,0, 0, 0, 0,0,0,0, 2];

    %  N = N'
%  M = M'
parse_gklf(N,M)
[V,info]=mpdefsub(M,N);
V
S = info.S;
M*V*S-N*V
eig(S)
[r,nu_r,nu_l] = nrankp(N,M)
