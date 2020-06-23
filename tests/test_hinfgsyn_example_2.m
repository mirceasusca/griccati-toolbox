%% Analytical solution
% Lecture_15_OCLS.pdf
% C:\Users\mircea\Desktop\Facultate\Master\Disertatie\bibliografie\robust\course
%
clf

P=tf([1 -10],[ 1 11 10]); [Ap,Bp,Cp,Dp]=ssdata(P);
np=max(size(Ap));
% filters for performance index
We=tf(1,[1 0.01]); [Ae,Be,Ce,De]=ssdata(We);
ne=max(size(Ae));
Wu=tf([1 2],[1 10]); [Au,Bu,Cu,Du]=ssdata(Wu);
nu=max(size(Au));
% augmented plant
A = [Ap, zeros(np,ne), zeros(np,nu);...
-Be*Cp, Ae, zeros(ne,nu);...
zeros(nu,np), zeros(nu,ne), Au];
B2 =[Bp; zeros(ne,1); Bu];
B1 =[Bp, zeros(np,1); zeros(ne,1), Be; zeros(nu,1), zeros(nu,1)];
C1 =[zeros(1,np), Ce, zeros(1,nu);...
zeros(1,np), zeros(1,ne), Cu];
D12=[ zeros(1,nu); Du];
C2 =[-Cp, zeros(1,ne), zeros(1,nu)];
D21=[0, 1];

Pg = ss(A,[B1,B2],[C1;C2],[zeros(size(D12,1),size(D21,2)),...
    D12;D21,zeros(size(D21,1),size(D12,2))]);
% Pg2 = augw(P,We,Wu,[]);

gamTry = 10.3;
[KM,CLM,GAMM] = hinfsyn(Pg,1,1,gamTry);
[K,CL,GAM,INFO] = hinfgsyn(Pg,1,1,gamTry);

GAMM
GAM

bode(K,KM), shg

%%
P = ss(P);
P.InputName = 'u';
P.OutputName = 'y1';
%
K.InputName = 'y2';
K.OutputName = 'u';
S = sumblk('y2 = r - y1');
CLP = connect(P,K,S,{'r'},{'y1'});

figure
subplot(121)
step(P,CLP)

subplot(122)
bode(P,CLP), shg