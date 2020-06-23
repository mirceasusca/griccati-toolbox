%% isstabilizable
rng(8);
n = 5;
m = 2;
p = 3;
%
A = rand(n);
B = rand(n,m);
C = rand(p,n);
D = rand(p,m);

% rank(ctrb(A,B))

assert(has_pos_sd_stab_sol(A,B,C,D) == 1)

%% works
A = 1;
B = 1;
C = 0;
D = 1;

assert(has_pos_sd_stab_sol(A,B,C,D) == 1)

%%
A = 1;
B = 1;
C = 1;
D = 0;  % fails at condition c)

assert(has_pos_sd_stab_sol(A,B,C,D) == 0)

%% works
A = -1;
B = 0;
C = 1;
D = 1;  % uncontrollable, but it has no unstable modes

assert(has_pos_sd_stab_sol(A,B,C,D) == 1)

%% works
A = [1,0;0,-1];
B = [0;1]; % unstabilizable
C = [1,1];
D = 0.00001; % invertible

assert(has_pos_sd_stab_sol(A,B,C,D) == 0)

%% works
A = [1,0;0,-1];
B = [1;0]; % stabilizable
C = [1,1];
D = 0.00001; % invertible

assert(has_pos_sd_stab_sol(A,B,C,D) == 1)

%% isstabilizable
rng(8);
n = 5;
m = 2;
p = 3;
%
A = rand(n);
B = rand(n,m);
C = rand(p,n);
D = rand(p,m);

discr = true;
% rank(ctrb(A,B))

assert(has_pos_sd_stab_sol(A,B,C,D,discr) == 1)

%% 
A = [0.435,0;0,-10.1];
B = [1;0]; % unstabilizable
C = [1,1];
D = 0.00001; % invertible

% sys = ss(A,B,C,D,1);
% t = 0:20;
% u = ones(1,length(t));
% lsim(sys,u,t,[0,0.1])

discr = true;
assert(has_pos_sd_stab_sol(A,B,C,D,discr) == 0)

%% 
A = [0.435,0;0,-10.1];
B = [0;1]; % uncontrollable, but stabilizable
C = [1,1];
D = 0.00001; % invertible

sys = ss(A,B,C,D,1);
% t = 0:20;
% u = ones(1,length(t));
% lsim(sys,u,t,[0,0.1])

discr = true;
assert(has_pos_sd_stab_sol(A,B,C,D,discr) == 1)