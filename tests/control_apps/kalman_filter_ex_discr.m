A = [1.1269   -0.4940    0.1129,
     1.0000         0         0,
          0    1.0000         0];

B = [-0.3832
      0.5919
      0.5191];

C = [1 0 0];

D = 0;

Plant = ss(A,[B B],C,0,-1,'inputname',{'u' 'w'},'outputname','y');
Q = 2; % A number greater than zero
R = 1; % A number greater than zero

[kalmf,L,X1,M,Z] = kalman(Plant,Q,R);
[L2,X2,info]=dlqe_df(A,B,C,0,Q,R)

%%
% First, build a complete plant model with u,w,v as inputs and
% y and yv as outputs:
a = A;
b = [B B 0*B];
c = [C;C];
d = [0 0 0;0 0 1];
P = ss(a,b,c,d,-1,'inputname',{'u' 'w' 'v'},'outputname',{'y' 'yv'});

sys = parallel(P,kalmf,1,1,[],[]);

SimModel = feedback(sys,1,4,2,1);
SimModel = SimModel([1 3],[1 2 3]);     % Delete yv form I/O

t = (0:300-1)';
u = sin(t/5);

rng(2);
w = sqrt(Q)*randn(length(t),1);
v = sqrt(R)*randn(length(t),1);

out = lsim(SimModel,[w,v,u]);

y = out(:,1);   % true response
ye = out(:,2);  % filtered response
ym = y + v;     % measured response

clf
subplot(211), plot(t,y,'b',t,ye,'r--',t,ym,'k-.'), grid minor
xlabel('No. of samples'), ylabel('Output') 
legend('y','ye','ym');
title('Kalman filter response')
subplot(212), plot(t,y-ym,'g',t,y-ye,'r--','linewidth',1.5), grid minor
xlabel('No. of samples'), ylabel('Error')
legend('y-ym','y-ye');

MeasErr = y-ym;
MeasErrCov = sum(MeasErr.*MeasErr)/length(MeasErr)
EstErr = y-ye;
EstErrCov = sum(EstErr.*EstErr)/length(EstErr)