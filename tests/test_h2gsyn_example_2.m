%% Analytical solution
% Lecture_14_OCLS.pdf
% C:\Users\mircea\Desktop\Facultate\Master\Disertatie\bibliografie\robust\course
%

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

% Find F
QQ = C1'*C1-C1'*D12*inv(D12'*D12)*D12'*C1;
AA = A-B2*inv(D12'*D12)*D12'*C1;
RR = B2*inv(D12'*D12)*B2';
Pu = are(AA,RR,QQ);
F = -inv(D12'*D12)*D12'*C1 - inv(D12'*D12)*B2'*Pu;
% Find L,
% C1->B1', D12->D21', A->A', B2->C2', F->L'
QQ = B1*B1'-B1*D21'*inv(D21*D21')*D21*B1';
AA = A'-C2'*inv(D21*D21')*D21*B1';
RR = C2'*inv(D21*D21')*C2;
Py = are(AA,RR,QQ);
L = ( -inv(D21*D21')*D21*B1' - inv(D21*D21')*C2*Py )';

% Compute the optimal controller
Ac=A+B2*F+L*C2; Bc=-L; Cc=F; Dc=0;
% K = minreal( ss(Ac,Bc,Cc,Dc) );
Ka = ss(Ac,Bc,Cc,Dc);

%%
Pg = ss(A,[B1,B2],[C1;C2],[zeros(size(D12,1),size(D21,2)),...
    D12;D21,zeros(size(D21,1),size(D12,2))]);
% Pg2 = augw(P,We,Wu,[]);

[KM,CLM,GAMM] = h2syn(Pg,1,1);
[K,CL,GAM,INFO] = h2gsyn(Pg,1,1);

GAMM
GAM

bode(Ka,KM,K), shg

%%
P = ss(P);
P.InputName = 'u';
P.OutputName = 'y1';
%
K.InputName = 'y2';
K.OutputName = 'u';
S = sumblk('y2 = r - y1');
CLP = connect(P,K,S,{'r'},{'y1'});
step(P,CLP)
% bode(P,CLP), shg