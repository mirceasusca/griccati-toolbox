rng(2); % rng=1 and rss=5 generates system with integrator; cont
sys = rss(5);
A = sys.A;
B = sys.B;
C = sys.C;
D = sys.D;
% 1.954039753986615
gamma = 1.95404; % modify it near hinfn = 1.954
% gamma = 1.95403;

hinfn = hinfnorm(sys,1e-7)
bounded_real_lemma(A,B,C,D,gamma)

%%
rng(1); % rng=1 and rss=5 generates system with integrator; cont
sys = rss(5);
A = sys.A;
B = sys.B;
C = sys.C;
D = sys.D;
gamma = 3.2;
% gamma = 2;

hinfn = hinfnorm(sys,1e-2)
bounded_real_lemma(A,B,C,D,gamma)

%% continuous
rng(8)
n = 5;
m = 3;
p = 2;
A = rand(n);
B = rand(n,m);
C = rand(p,n);
D = rand(p,m);
gamma = 2;
sys = ss(A,B,C,D);
hinfn = hinfnorm(sys,1e-2)
bounded_real_lemma(A,B,C,D,gamma)

%% discrete system
rng(2); % rng=1 and rss=5 generates system with integrator
sys = drss(5);
A = sys.A;
B = sys.B;
C = sys.C;
D = sys.D;
% gamma = 2.05; % check computations and tolerances
gamma = 2.086693;

hinfn = hinfnorm(sys,1e-8)
discr = true;
bounded_real_lemma(A,B,C,D,gamma,discr)

%%
rng(1); % rng=1 and rss=5 generates system with integrator; cont
sys = drss(9);
A = sys.A;
B = sys.B;
C = sys.C;
D = sys.D;
gamma = 5.2;

discrete = true;
hinfn = hinfnorm(sys,1e-2)
bounded_real_lemma(A,B,C,D,gamma,discr)
