%% Nominal H-Infinity Design
% Use hinfsyn to compute an Hinf controller for each value of the blending
% vector beta
ncont = 1; % one control signal, u
nmeas = 2; % two measurement signals, sd and ab
Kinf = ss(zeros(ncont,nmeas,3));
gammainf = zeros(3,1);
for i=1:3
%    [Kinf(:,:,i),~,gammainf(i)] = hinfsyn(qcaric(:,:,i),nmeas,ncont);
   [Kinf(:,:,i),~,gammainf(i)] = hinfgsyn(qcaric(:,:,i),nmeas,ncont,1);
end
%
display(gammainf,'gamma Hinf');
%
% Construct the corresponding closed-loop models and compare the gains from
% road disturbance to xb, sd, ab for the passive and active suspensions.
% Observe that all three controllers reduce suspension deflection and body
% acceleration below the rattlespace frequency (23 rad/s)
Kinf.u = {'sd','ab'};  Kinf.y = 'u';
CLinf = connect(qcar,Act.Nominal,Kinf,'r',{'xb';'sd';'ab'});
%
figure;
bodemag(qcar(:,'r'),'b', CLinf(:,:,1),'r-.', ...
   CLinf(:,:,2),'m-.', CLinf(:,:,3),'k-.',{1,140}), grid
set(findall(gcf,'type','line'),'linewidth',1.25)
legend('Open-loop','Comfort','Balanced','Handling','location','SouthEast')
title('H_{\infty}: Body travel, suspension deflection, and body acceleration due to road')

%% Time-Domain Evaluation
% To further evaluate the three designs, perform time-domain simulations
% using a road disturbance signal r(t) representing a road bump of height 5
% cm.
% Road disturbance
t = 0:0.0025:1;
roaddist = zeros(size(t));
roaddist(1:101) = 0.025*(1-cos(8*pi*t(1:101)));
%
% Closed-loop model
SIMKinf = connect(qcar,Act.Nominal,Kinf,'r',{'xb';'sd';'ab';'fs'});

% Simulate
p1 = lsim(qcar(:,1),roaddist,t);
y1 = lsim(SIMKinf(1:4,1,1),roaddist,t);
y2 = lsim(SIMKinf(1:4,1,2),roaddist,t);
y3 = lsim(SIMKinf(1:4,1,3),roaddist,t);

% Plot results
figure;
subplot(211)
plot(t,roaddist,'g',t,p1(:,1),'b',t,y1(:,1),'r.',t,y2(:,1),'m.',t,y3(:,1),'k.','linewidth',1.5); grid;
title('H_{\infty}: Body travel'), ylabel('x_b (m)')
legend('Road','Open-loop','Comfort','Balanced','Handling');
subplot(212)
plot(t,roaddist,'g',t,p1(:,3),'b',t,y1(:,3),'r.',t,y2(:,3),'m.',t,y3(:,3),'k.','linewidth',1.5); grid;
title('H_{\infty}: Body acceleration'), xlabel('Time (s)'), ylabel('a_b (m/s^2)')
%
% Observe that the body acceleration is smallest for the controller
% emphasizing passenger comfort and largest for the controller emphasizing
% suspension deflection. The "balanced" design achieves a good compromise
% between body acceleration and suspension deflection.

%%
figure;
subplot(211)
plot(t,roaddist,'g',t,p1(:,2),'b',t,y1(:,2),'r.',t,y2(:,2),'m.',t,y3(:,2),'k.','linewidth',1.5); grid;
title('H_{\infty}: Suspension deflection'), xlabel('Time (s)'), ylabel('s_d (m)')
subplot(212)
plot(t,roaddist,'g',t,zeros(size(t)),'b',t,y1(:,4),'r.',t,y2(:,4),'m.',t,y3(:,4),'k.','linewidth',1.5); grid;
title('H_{\infty}: Control force'), xlabel('Time (s)'), ylabel('f_s (kN)')
legend('Road Disturbance','Open-loop','Comfort','Balanced','Handling','location','SouthEast')

%% Simulate the response to a road bump for 100 actuator models
% randomly selected from the uncertain model Act
rng('default'), nsamp = 100;

% Uncertain closed-loop model with balanced H-infinity controller
CLU = connect(qcar,Act,Kinf(:,:,2),'r',{'xb','sd','ab'});
figure;
lsim(usample(CLU,nsamp),'b',CLU.Nominal,'r',roaddist,t); grid;
title('H_{\infty}: Nominal "balanced" design')
legend('Perturbed','Nominal','location','SouthEast')
