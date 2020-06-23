%% isstabilizable - CTS
rng(8);
n = 5;
m = 3;
p = 2;
%
A = rand(n);
B = rand(n,m);
C = rand(p,n);
D = rand(p,m);

rank(ctrb(A,B))

isstabilizable(A,B,C,D)
isdetectable(A,B,C,D)

%% isstabilizable - DTS
rng(8);
n = 5;
m = 3;
p = 2;
%
A = rand(n);
B = rand(n,m);
C = rand(p,n);
D = rand(p,m);

rank(ctrb(A,B))

discr = true;
isstabilizable(A,B,C,D,discr)
isdetectable(A,B,C,D,discr)

%%
A = 1;
B = 1;
C = 0;
D = 0;

iss = isstabilizable(A,B,C,D)
isd = isdetectable(A,B,C,D)

%%
A = 1;
B = 0;
C = 1;
D = 0;

iss = isstabilizable(A,B,C,D)
isd = isdetectable(A,B,C,D)
