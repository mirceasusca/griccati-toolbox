% continuous system, stable feedback
n = 14;
m = 2;

M = rand(n);
N = rand(n);
B = rand(n,m);

opt.stable = true;
opt.discr = false;
sdeg = -0.35;
F = gsfstab(N,M,B,[],sdeg,opt)
real(eig(N+B*F,M))

%%
% continuous system, stable feedback
n = 14;
m = 2;

M = rand(n);
N = rand(n);
B = rand(n,m);

sdeg = -0.35;
F = gsfstab(N,M,B,[],sdeg)
real(eig(N+B*F,M))

%%
% continuous system, stable feedback
n = 14;
m = 2;

M = rand(n);
N = rand(n);
B = rand(n,m);

F = gsfstab(N,M,B)
real(eig(N+B*F,M))

%% continuous system, antistable feedback
n = 14;
m = 2;

M = rand(n);
N = rand(n);
B = rand(n,m);

opt.stable = false;
opt.discr = false;
sdeg = 0.35;
F = gsfstab(N,M,B,[],sdeg,opt)
real(eig(N+B*F,M))

%% discrete system, stable feedback
n = 14;
m = 2;

M = rand(n);
N = rand(n);
B = rand(n,m);

opt.stable = true;
opt.discr = true;
sdeg = 0.45;
F = gsfstab(N,M,B,[],sdeg,opt)
abs(eig(N+B*F,M))

%% discrete system, stable feedback
n = 14;
m = 2;

M = rand(n);
N = rand(n);
B = rand(n,m);

sdeg = 0.95;
F = gsfstab(N,M,B,[],sdeg)
abs(eig(N+B*F,M))

%% discrete system, antistable feedback
n = 14;
m = 2;

M = rand(n);
N = rand(n);
B = rand(n,m);

opt.stable = false;
opt.discr = true;
sdeg = 1.45;
F = gsfstab(N,M,B,[],sdeg,opt)
abs(eig(N+B*F,M))
