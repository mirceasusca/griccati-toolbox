%% Nominal H2 Design
% Use h2syn to compute an H2 controller for each value of the blending
% vector beta
ncont = 1; % one control signal, u
nmeas = 2; % two measurement signals, sd and ab
K2 = ss(zeros(ncont,nmeas,3));
gamma2 = zeros(3,1);
for i=1:3
%    [K2(:,:,i),~,gamma2(i)] = h2syn(ss(qcaric(:,:,i)),nmeas,ncont);
   [K2(:,:,i),~,gamma2(i)] = h2gsyn(ss(qcaric(:,:,i)),nmeas,ncont);
end
%
display(gamma2,'gamma H2');
%
% Construct the corresponding closed-loop models and compare the gains from
% road disturbance to xb, sd, ab for the passive and active suspensions.
% Observe that all three controllers reduce suspension deflection and body
% acceleration below the rattlespace frequency (23 rad/s)
K2.u = {'sd','ab'};  K2.y = 'u';
CL2 = connect(qcar,Act.Nominal,K2,'r',{'xb';'sd';'ab'});
%
figure;
bodemag(qcar(:,'r'),'b', CL2(:,:,1),'r-.', ...
   CL2(:,:,2),'m-.', CL2(:,:,3),'k-.',{1,140}), grid
set(findall(gcf,'type','line'),'linewidth',1.25)
legend('Open-loop','Comfort','Balanced','Handling','location','SouthEast')
title('H_2: Body travel, suspension deflection, and body acceleration due to road')

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
SIMK2 = connect(qcar,Act.Nominal,K2,'r',{'xb';'sd';'ab';'fs'});

% Simulate
p1 = lsim(qcar(:,1),roaddist,t);
y1 = lsim(SIMK2(1:4,1,1),roaddist,t);
y2 = lsim(SIMK2(1:4,1,2),roaddist,t);
y3 = lsim(SIMK2(1:4,1,3),roaddist,t);

% Plot results
figure;
subplot(211)
plot(t,roaddist,'g',t,p1(:,1),'b',t,y1(:,1),'r.',t,y2(:,1),'m.',t,y3(:,1),'k.','linewidth',1.5); grid;
title('H_2: Body travel'), ylabel('x_b (m)')
legend('Road','Open-loop','Comfort','Balanced','Handling');
subplot(212)
plot(t,roaddist,'g',t,p1(:,3),'b',t,y1(:,3),'r.',t,y2(:,3),'m.',t,y3(:,3),'k.','linewidth',1.5); grid;
title('H_2: Body acceleration'), xlabel('Time (s)'), ylabel('a_b (m/s^2)')
%
% Observe that the body acceleration is smallest for the controller
% emphasizing passenger comfort and largest for the controller emphasizing
% suspension deflection. The "balanced" design achieves a good compromise
% between body acceleration and suspension deflection.

%%
figure;
subplot(211)
plot(t,roaddist,'g',t,p1(:,2),'b',t,y1(:,2),'r.',t,y2(:,2),'m.',t,y3(:,2),'k.','linewidth',1.5); grid;
title('H_2: Suspension deflection'), xlabel('Time (s)'), ylabel('s_d (m)')
subplot(212)
plot(t,roaddist,'g',t,zeros(size(t)),'b',t,y1(:,4),'r.',t,y2(:,4),'m.',t,y3(:,4),'k.','linewidth',1.5); grid;
title('H_2: Control force'), xlabel('Time (s)'), ylabel('f_s (kN)')
legend('Road Disturbance','Open-loop','Comfort','Balanced','Handling','location','SouthEast')

%% Simulate the response to a road bump for 100 actuator models
% randomly selected from the uncertain model Act
rng('default'), nsamp = 100; 
%
% Uncertain closed-loop model with balanced H2 controller
CLU = connect(qcar,Act,K2(:,:,2),'r',{'xb','sd','ab'});
figure;
lsim(usample(CLU,nsamp),'b',CLU.Nominal,'r',roaddist,t); grid;
title('H_2: Nominal "balanced" design')
legend('Perturbed','Nominal','location','SouthEast')
