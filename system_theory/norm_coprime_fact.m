function [M,N,U,W,Mt,Nt,Ut,Wt,info] = norm_coprime_fact(A,B,C,Ts)
%NORM_COPRIME_FACT Given the system G = ss(A,B,C,0,Ts), the pair (A,B) is
% stabilizable and the pair (C,A) is detectable, the function returns the
% realizations in terms of normalized coprime factorizations of a 
% continuous-time or discrete-time system:
%  G = N*M^-1, with N'*N + M'*M = I, and
%  G = Mt^-1*Nt, with Nt*Nt' + Mt*Mt' = I.
%
% [M,N,U,W,Mt,Nt,Ut,Wt,info] = norm_coprime_fact(A,B,C), default Ts=0.
% [M,N,U,W,Mt,Nt,Ut,Wt,info] = norm_coprime_fact(A,B,C,Ts)
%

if nargin == 3
    Ts = 0;
end

assert(size(A,1) == size(B,1),'Incompatible system dimensions');
assert(size(A,2) == size(C,2),'Incompatible system dimensions');

[n,m] = size(B);
[p,~] = size(C);

sigma_x = create_popov_triplet(A,B,C'*C,zeros(n,m),eye(m),Ts ~= 0);
sigma_y = create_popov_triplet(A',C',B*B',zeros(n,p),eye(p),Ts ~= 0);

[F,X,info_x]=gricsv(sigma_x);
[K,Y,info_y]=gricsv(sigma_y);
K = K';

info.info_x = info_x;
info.info_y = info_y;

AF = A + B*F;
AK = A + K*C;

if Ts == 0
    M = ss(AF,B,F,eye(m));
    N = ss(AF,B,C,zeros(p,m));
    U = ss(AK,B,-F,eye(m));
    W = ss(AK,K,F,zeros(m,p));

    Mt = ss(AK,K,C,eye(p));
    Nt = ss(AK,B,C,zeros(p,m));
    Ut = ss(AF,-K,C,eye(p));
    Wt = ss(AF,K,F,zeros(m,p));
else
    V = chol(eye(m) + B'*X*B);
    Vt = chol(eye(p) + C*Y*C'); Vt = Vt';
    
    M = ss(AF,B/V,F,inv(V),Ts);
    N = ss(AF,B/V,C,zeros(p,m),Ts);
    U = ss(AK,B,-V*F,V,Ts);
    W = ss(AK,K,V*F,zeros(m,p),Ts);
    Mt = ss(AK,K,Vt\C,inv(Vt),Ts);
    Nt = ss(AK,B,Vt\C,zeros(p,m),Ts);
    Ut = ss(AK,-K*Vt,C,Vt,Ts);
    Wt = ss(AF,K*Vt,F,zeros(m,p),Ts);
end

% for debug
% G = ss(A,B,C,0,Ts);
% G1 = series(inv_sys(M),N);
% G2 = series(Nt,inv_sys(Mt));
% bode(G,G1,G2);

end