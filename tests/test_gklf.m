%% Example from the authors
A = [0 1 0; 0 0 2];
B = [0 0 0; 0 0 3];
[At,Bt,info,Q,Z] = gklf(A,B,0)
% Q*At*Z'-A
% Q*Bt*Z'-B
parse_gklf(A,B)
% norm(P'*A*Q-S)
% norm(P'*B*Q-T)

%% Regular pencil with only finite generalized eigenvalues
rng(3);
A = rand(3);
B = rand(3);
% [S,T,P,Q,kstr]= guptri(A,B);
% lambda1 = eig(A,B)
% lambda2 = diag(S)./diag(T)
[At,Bt,info,Q,Z] = gklf(A,B,0);
parse_gklf(A,B)
% Q*At*Z'-A
% Q*Bt*Z'-B
% 
%% Regular pencil with only infinite generalized eigenvalues
% rng(3);
A = diag([1,2,3]);
B = zeros(3);
B(1,2) = 1;

[At,Bt,info,Q,Z] = gklf(A,B,0);
% Q*At*Z'-A
% Q*Bt*Z'-B
parse_gklf(A,B)

%% Regular pencil with only finite and infinite generalized eigenvalues
% rng(3);
A = diag([1,2,3]);
B = [0,0,0;0,1,0;0,0,1];

rng(5)
Q = qr(rand(3));
A = Q*A*Q';
B = Q*B*Q';

[At,Bt,info,Q,Z] = gklf(A,B,0);
parse_gklf(A,B)
% Q*At*Z'-A
% Q*Bt*Z'-B

%% Regular pencil with only zero eigenvalues
% rng(3);
A = zeros(3);
B = diag([1,2,3]);
% [S,T,P,Q,kstr]= guptri(A,B);
% 
[At,Bt,info,Q,Z] = gklf(A,B,0);
parse_gklf(A,B)
% Q*At*Z'-A
% Q*Bt*Z'-B

%% Regular pencil with finite and infinite generalized eigenvalues
rng(3);
[Q,~] = qr(rand(3));
A = Q*diag([1,2,3])*Q';
B = [0, 0, 0; 0, 1, 0; 0, 0, 1];
[At,Bt,info,Q,Z] = gklf(A,B,0);
parse_gklf(A,B)
% parse_kstr(kstr)
% isregular(A,B)

%% Singular pencil with right Kronecker indices
A = diag(ones(3,1),1); A(4,:) = [];
B = eye(3,4);
% [S,T,P,Q,kstr]= guptri(A,B);

[At,Bt,info,Q,Z] = gklf(A,B,0)
% Q*At*Z'-A
% Q*Bt*Z'-B

%% Singular pencil with left Kronecker indices
A = diag(ones(3,1),1); A(4,:) = []; A = A';
B = eye(3,4); B = B';

[At,Bt,info,Q,Z] = gklf(A,B,0);
parse_gklf(A,B)
% Q*At*Z'-A
% Q*Bt*Z'-B

%% Example with all 5 cases: right indices, zero, finite, 
% infinite, left indices
A = [
    22  34  31   31  17 
    45  45  42   19  29  
    39  47  49   26  34 
    27  31  26   21  15 
    38  44  44   24  30 
    ];

B = [
    13  26  25  17  24 
    31  46  40  26  37 
    26  40  19  25  25 
    16  25  27  14  23 
    24  35  18  21  22  
    ];

% [S,T,P,Q,kstr]= guptri(A,B)
[At,Bt,info,Q,Z] = gklf(A,B,0);
parse_gklf(A,B)
% Q*At*Z'-A
% Q*Bt*Z'-B

%% Example with multiple singularities
% infinite, left indices
A = [
    22  34  31   31  17 34
    22  34  31   31  17 34
    45  45  42   19  29 45 
    39  47  49   26  34 47
    27  31  26   21  15 31
    38  44  44   24  30 44
    ];

B = [
    22  34  31   31  17 34
    13  26  25  17  24 26
    31  46  40  26  37 46
    26  40  19  25  25 40
    16  25  27  14  23 25
    24  35  18  21  22 35 
    ];

% [S,T,P,Q,kstr]= guptri(A,B);
% parse_kstr(kstr)
% isregular(A,B)

[At,Bt,info,Q,Z] = gklf(A,B,0);
parse_gklf(A,B)
% Q*At*Z'-A
% Q*Bt*Z'-B
% info

%% Example with two finite eigenvalues and a left Kronecker index
A = zeros(4,3);
B = zeros(4,3);
A(1,1) = 1;
A(2,2) = 1;
A(4,3) = 1;
B(1,1) = 1;
B(2,2) = 2;
B(3,3) = 1;

[At,Bt,info,Q,Z] = gklf(A,B,0);
parse_gklf(A,B);
% Q*At*Z'-A
% Q*Bt*Z'-B

%%
A = zeros(4);
B = zeros(4);
A(1,1) = 1;
A(2,2) = 1;
A(4,4) = 1;
B(1,2) = 1;
B(2,2) = 2;
B(3,4) = 1;

% [S,T,P,Q,kstr]= guptri(A,B);
% parse_kstr(kstr)
[At,Bt,info,Q,Z] = gklf(A,B,0);
parse_gklf(A,B);
% Q*At*Z'-A
% Q*Bt*Z'-B

%% Example with multiple singularities
A = zeros(8,12);
B = zeros(8,12);
B(1,1) = 1;
B(2,2) = 1;
B(3,3) = 1;
B(4,4) = 1;
B(5,7) = 1;
B(6,9) = 1;
B(7,10) = 1;
B(8,11) = 1;
A(2,3) = 1;
A(3,4) = 1;
A(5,8) = 1;
A(6,10) = 1;
A(7,11) = 1;
A(8,12) = 1;

[At,Bt,info,Q,Z] = gklf(A,B);
parse_gklf(A,B);
% Q*At*Z'-A
% Q*Bt*Z'-B

%% Regular pencil with finite and infinite generalized eigenvalues
rng(3);
N = 8;
[Q1,~] = qr(rand(8));
[Q2,~] = qr(rand(8));
A = Q1*diag([1:8])*Q1';
B = zeros(8); B(4:8,4:8) = eye(5); B = Q2'*B*Q2;
%
[At,Bt,info,Q,Z] = gklf(A,B);
parse_gklf(A,B);
% [S,T,P,Q,kstr]= guptri(A,B);
% lambda1 = eig(A,B)
% lambda2 = diag(S)./diag(T); lambda2(~isfinite(lambda2)) = Inf
% parse_kstr(kstr)
% isregular(A,B)

%% Regular pencil with finite and infinite generalized eigenvalues
rng(3);
N = 14;
[Q1,~] = qr(rand(N));
[Q2,~] = qr(rand(N));
A = diag([1:N]); A(4:5,4:5) = [1,-2;2,1]; A = Q1*A*Q1';
B = zeros(N); B(1:4,1:4) = eye(4); B = Q2'*B*Q2;

[At,Bt,info,Q,Z] = gklf(A,B);
parse_gklf(A,B);
% lambda1 = eig(A,B)
% lambda2 = diag(S)./diag(T); lambda2(~isfinite(lambda2)) = Inf

%% Regular pencil with finite and infinite generalized eigenvalues
rng(3);
N = 14;
[Q1,~] = qr(rand(N));
[Q2,~] = qr(rand(N));
A = diag([1:round(N/2),-(round(N/2)+1:N)]); A(4:5,4:5) = [-1,-2;2,-1];%  A = Q1*A*Q1';
B = zeros(N); B(1:5,1:5) = eye(5); B(8:13,8:13) = eye(13-8+1); % B = Q2'*B*Q2; 
% [S,T,P,Q,kstr]= guptri(A,B);
% parse_kstr(kstr)
% isregular(A,B)
[At,Bt,info,Q,Z] = gklf(A,B,0);
parse_gklf(A,B,0);
% Q*At*Z'-A
% Q*Bt*Z'-B
% lambda1 = eig(A,B)
% lambda2 = diag(S)./diag(T); lambda2(~isfinite(lambda2)) = Inf
% [SS,TT,QQ,ZZ] = ordqz(S,T,P,Q,'lhp')

%%
A = diag(ones(1,8)); A(1,1) = 0; A(2,2) = 0;
B = diag([1,2,3,4,5,0,0,0]);
% [S,T,P,Q,kstr]= guptri(A,B);

[At,Bt,info,Q,Z] = gklf(A,B,0);
parse_gklf(A,B)

%% L1,2(lambda) right Kron.
A = [0,1];
B = [1,0];
% [S,T,P,Q,kstr]= guptri(A,B);

[At,Bt,info,Q,Z] = gklf(A,B,0)
parse_gklf(A,B)

%% L2,1(lambda) left Kron.
A = [0;1];
B = [1;0];
% [S,T,P,Q,kstr]= guptri(A,B);

[At,Bt,info,Q,Z] = gklf(B,A,0)
parse_gklf(A,B,0)