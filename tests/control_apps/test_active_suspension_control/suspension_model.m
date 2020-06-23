%%
clear all
clc
close all

%% Physical parameters
mb = 300;    % kg
mw = 60;     % kg
bs = 1000;   % N/m/s
ks = 16000 ; % N/m
kt = 190000; % N/m

% State matrices
A = [ 0 1 0 0; [-ks -bs ks bs]/mb ; ...
      0 0 0 1; [ks bs -ks-kt -bs]/mw];
B = [ 0 0; 0 1e3/mb ; 0 0 ; [kt -1e3]/mw];
C = [1 0 0 0; 1 0 -1 0; A(2,:)];
D = [0 0; 0 0; B(2,:)];

qcar = ss(A,B,C,D);
qcar.StateName = {'body travel (m)';'body vel (m/s)';...
          'wheel travel (m)';'wheel vel (m/s)'};
qcar.InputName = {'r';'fs'};
qcar.OutputName = {'xb';'sd';'ab'};

%% Invariant zeros of the system
tzero(qcar({'xb','ab'},'fs'))
% imaginary-axis zero with natural frequency 56.27 rad/s -- tire-hop freq.

%% Zero of the system
zero(qcar('sd','fs'))
% imaginary-axis zero with natural frequency 22.97 rad/s -- rattlespace fr.

%%
% Road disturbances influence the motion of the car and suspension.
% Passenger comfort is associated with small body acceleration. The
% allowable suspension travel is constrained by limits of the actuator
% displacement.
figure;
bodemag(qcar({'ab','sd'},'r'),'b',qcar({'ab','sd'},'fs'),'r',{1,100});
set(findall(gcf,'type','line'),'linewidth',1.5);
legend('Road disturbance (r)', 'Actuator force (fs)', 'location', 'southwest');
title(['Body acceleration (ab) and suspension travel (sd) gains']); grid;

%% Step response
figure
subplot(121);
step(qcar({'xb','ab','sd'},'r'),'b',3);
set(findall(gcf,'type','line'),'linewidth',1.25)
title(['Step response from road disturbance (r)']); grid;
%
% figure
subplot(122)
step(qcar({'xb','ab','sd'},'fs'),'r',3);
set(findall(gcf,'type','line'),'linewidth',1.25)
% legend('Road disturbance input (r)', 'Actuator force input (fs)', 'location', 'southwest');
title('Step response from actuator force (fs)'); grid;

%% Uncertain Actuator Model (hydraulic actuator)
% family of actuator odels to account for modeling errors and variability
% in the acttuator and quarter-car models
% Nominal model + frequency-dependent amount of uncertainty (modulated)
ActNom = tf(1,[1/60,1]);
%
% 1st ord sys with given DC, CROSSW, HF
% G(j*0)=DC, |G(j*CROSSW)|=1, G(j*infty)=HF
Wunc = makeweight(0.40,15,3);
% Create a block unc which represents arbitrary LTI dynamics with frequency
% response gain no larger than one.
unc = ultidyn('unc',[1,1],'SampleStateDim',5);
Act = ActNom*(1+Wunc*unc);
Act.InputName='u';
Act.OutputName='fs';
%
% rng('default');
figure;
bode(Act,'b',Act.NominalValue,'r+',logspace(-1,3,120)); grid;
% set(findall(gcf,'type','line'),'linewidth',1.25)

%% Design Setup
% measurement y1 of the suspension travel sd
% measurement y2 of the body acceleration ab
% compute the control signal u driving the hydraulic actuator
% disturbance 1: road disturbance r, modeled as a normalized signal d1
% shaped by a weighting function Wroad
% sensor noise on both measurements: normalized signals d2 and d3 shaped by
% weighting functions Wd2 and Wd3.
%
% Objectives: passenger comfort (ab) and road handling (sd);
% Disturbance rejection goal: minimize the impact of the disturbances d1,
% d2, d3 on a weighted combination of control effort u, suspension travel
% sd, and body acceleration ab.
%
% When using Hinf norm (peak gain) to measure "impact", this amounts to
% designing a controller that minimizez Hinf norm from disturbance inputs
% d1, d2, d3 to error signals e1, e2, e3.
%
Wact = 0.8*tf([1 50],[1 500]);  Wact.u = 'u';  Wact.y = 'e1';
Wroad = ss(0.07);  Wroad.u = 'd1';   Wroad.y = 'r';
Wd2 = ss(0.01);  Wd2.u = 'd2';   Wd2.y = 'Wd2';
Wd3 = ss(0.5);   Wd3.u = 'd3';   Wd3.y = 'Wd3';
%
% Specify closed-loop targets for the gain from road disturbance r to
% suspension deflection sd (handling) and body acceleration ab (comfort).
% Because of the actuator uncertainty and imaginary-axis zeros, only seek
% to attenuate disturbances below 10 rad/s.
HandlingTarget = 0.04 * tf([1/8 1],[1/80 1]);
ComfortTarget = 0.4 * tf([1/0.45 1],[1/150 1]);
Targets = [HandlingTarget ; ComfortTarget];
%
figure;
bodemag(qcar({'sd','ab'},'r')*Wroad,'b',Targets,'r--',{1,1000}), grid
set(findall(gcf,'type','line'),'linewidth',1.5);
title('Response to road disturbance')
legend('Open-loop','Closed-loop target')

%% Three design points
% The corresponding performance weights Wsd, Wab are the reciprocals of the
% comfort and handling targets. To investigate the trade-off between
% passenger comfort and road handling, construct three sets of weights
% (beta*Wsd, (1-beta)*Wab) corresponding to three trade-offs: comfort
% (beta=0.01), balanced (beta=0.5), handling (beta=0.99)
beta = reshape([0.01 0.5 0.99],[1 1 3]);
Wsd = beta / HandlingTarget;
Wsd.u = 'sd';  Wsd.y = 'e3';
Wab = (1-beta) / ComfortTarget;
Wab.u = 'ab';  Wab.y = 'e2';
% Construct a model qcaric of the Fig. 2 block diagram.
% qcaric is an array of three models, one for each design point beta. Also,
% it is an uncertain model since it contains the uncertain actuator model
% Act.
sdmeas  = sumblk('y1 = sd+Wd2');
abmeas = sumblk('y2 = ab+Wd3');
ICinputs = {'d1';'d2';'d3';'u'};
ICoutputs = {'e1';'e2';'e3';'y1';'y2'};
qcaric = connect(qcar(2:3,:),Act,Wroad,Wact,Wab,Wsd,Wd2,Wd3,...
                 sdmeas,abmeas,ICinputs,ICoutputs);
display(qcaric);
