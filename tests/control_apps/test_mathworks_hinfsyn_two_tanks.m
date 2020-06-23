close all, clc, clearvars

A1 = 0.0256;	% Area of tank 1 (hunits^2)
A2 = 0.0477;	% Area of tank 2 (hunits^2)
h2 = 0.241;	    % Height of tank 2, fixed by overflow (hunits)
fb = 3.28e-5;   % Bias stream flow (hunits^3/sec)
fs = 0.00028;	% Flow scaling (hunits^3/sec/funit)
th = 1.0;	    % Hot water supply temp (tunits)
tc = 0.0;	    % Cold water supply temp (tunits)
tb = tc;	    % Cold bias stream temp (tunits)
alpha = 4876;   % Constant for flow/height relation (hunits/funits)
beta = 0.59;    % Constant for flow/height relation (hunits)

h1ss = 0.75;                            % Water level for tank 1
t1ss = 0.75;                            % Temperature of tank 1
f1ss = (h1ss+beta)/alpha;               % Flow tank 1 -> tank 2
fss = [th,tc;1,1]\[t1ss*f1ss;f1ss];
fhss = fss(1);                          % Hot flow
fcss = fss(2);                          % Cold flow
t2ss = (f1ss*t1ss + fb*tb)/(f1ss + fb); % Temperature of tank 2

A = [ -1/(A1*alpha),          0;
      (beta*t1ss)/(A1*h1ss),  -(h1ss+beta)/(alpha*A1*h1ss)];

B = fs*[ 1/(A1*alpha),   1/(A1*alpha);
         th/A1,          tc/A1];

C = [ alpha,             0;
      -alpha*t1ss/h1ss,  1/h1ss];

D = zeros(2,2);
tank1nom = ss(A,B,C,D,'InputName',{'fh','fc'},'OutputName',{'h1','t1'});

figure
step(tank1nom), title('Step responses of Tank 1')

A = -(h1ss + beta + alpha*fb)/(A2*h2*alpha);
B = [ (t2ss+t1ss)/(alpha*A2*h2),  (h1ss + beta)/(alpha*A2*h2) ];
C = 1;
D = zeros(1,2);

tank2nom = ss(A,B,C,D,'InputName',{'h1','t1'},'OutputName','t2');

% figure
step(tank2nom), title('Step responses of Tank 2')

%% Actuator models
act_BW = 20;		% Actuator bandwidth (rad/sec)
actuator = [ tf(act_BW,[1 act_BW]); tf([act_BW 0],[1 act_BW]) ];
actuator.OutputName = {'Flow','Flow rate'};

bodemag(actuator)
title('Valve actuator dynamics')

hot_act = actuator;
set(hot_act,'InputName','fhc','OutputName',{'fh','fh_rate'});
cold_act =actuator;
set(cold_act,'InputName','fcc','OutputName',{'fc','fc_rate'});

% Anti-aliasing filters
fbw = 2.25;		% Anti-aliasing filter cut-off (Hz)
filter = mkfilter(fbw,4,'Butterw');
h1F = filter;
t1F = filter;

%% Uncertainty on Model Dynamics
Wh1 = 0.01+tf([0.5,0],[0.25,1]);
Wt1 = 0.1+tf([20*h1ss,0],[0.2,1]);
Wt2 = 0.1+tf([100,0],[1,21]);

clf
bodemag(Wh1,Wt1,Wt2), title('Relative bounds on modeling errors')
legend('h1 dynamics','t1 dynamics','t2 dynamics','Location','NorthWest')
t2F = filter;

% Normalized error dynamics
delta1 = ultidyn('delta1',[1 1]);
delta2 = ultidyn('delta2',[1 1]);
delta3 = ultidyn('delta3',[1 1]);

% Frequency-dependent variability in h1, t1, t2 dynamics
varh1 = 1+delta1*Wh1;
vart1 = 1+delta2*Wt1;
vart2 = 1+delta3*Wt2;

% Add variability to nominal models
tank1u = append(varh1,vart1)*tank1nom;
tank2u = vart2*tank2nom;

tank1and2u = [0 1; tank2u]*tank1u;

% figure;
step(tank1u,1000), title('Variability in responses due to modeling errors (Tank 1)')

%% Setting up a controller design
Wh1noise = zpk(0.01);  % h1 noise weight
Wt1noise = zpk(0.03);  % t1 noise weight
Wt2noise = zpk(0.03);  % t2 noise weight

Wt1perf = tf(100,[400,1]);	% t1 tracking error weight
Wt2perf = tf(50,[800,1]);	% t2 tracking error weight

clf
bodemag(Wt1perf,Wt2perf)
title('Frequency-dependent penalty on setpoint tracking errors')
legend('t1','t2')

Wt1cmd = zpk(0.1);               % t1 input command weight
Wtdiffcmd = zpk(0.01);           % t2 - t1  input command weight

Whact =  zpk(0.01);  % Hot actuator penalty
Wcact =  zpk(0.01);  % Cold actuator penalty

Whrate = zpk(50);    % Hot actuator rate penalty
Wcrate = zpk(50);    % Cold actuator rate penalty

%% Building a Weighted Open-Loop Model
inputs = {'t1cmd', 'tdiffcmd', 't1noise', 't2noise', 'fhc', 'fcc'};
outputs = {'y_Wt1perf', 'y_Wt2perf', 'y_Whact', 'y_Wcact', ...
             'y_Whrate', 'y_Wcrate', 'y_Wt1cmd', 'y_t1diffcmd', ...
                                           'y_t1Fn', 'y_t2Fn'};

hot_act.InputName = 'fhc'; hot_act.OutputName = {'fh' 'fh_rate'};
cold_act.InputName = 'fcc'; cold_act.OutputName = {'fc' 'fc_rate'};

tank1and2u.InputName = {'fh','fc'};
tank1and2u.OutputName = {'t1','t2'};

t1F.InputName = 't1'; t1F.OutputName = 'y_t1F';
t2F.InputName = 't2'; t2F.OutputName = 'y_t2F';

Wt1cmd.InputName = 't1cmd'; Wt1cmd.OutputName = 'y_Wt1cmd';
Wtdiffcmd.InputName = 'tdiffcmd'; Wtdiffcmd.OutputName = 'y_Wtdiffcmd';

Whact.InputName = 'fh'; Whact.OutputName = 'y_Whact';
Wcact.InputName = 'fc'; Wcact.OutputName = 'y_Wcact';

Whrate.InputName = 'fh_rate'; Whrate.OutputName = 'y_Whrate';
Wcrate.InputName = 'fc_rate'; Wcrate.OutputName = 'y_Wcrate';

Wt1perf.InputName = 'u_Wt1perf'; Wt1perf.OutputName = 'y_Wt1perf';
Wt2perf.InputName = 'u_Wt2perf'; Wt2perf.OutputName = 'y_Wt2perf';

Wt1noise.InputName = 't1noise'; Wt1noise.OutputName = 'y_Wt1noise';
Wt2noise.InputName = 't2noise'; Wt2noise.OutputName = 'y_Wt2noise';

sum1 = sumblk('y_t1diffcmd = y_Wt1cmd + y_Wtdiffcmd');
sum2 = sumblk('y_t1Fn = y_t1F + y_Wt1noise');
sum3 = sumblk('y_t2Fn = y_t2F + y_Wt2noise');
sum4 = sumblk('u_Wt1perf = y_Wt1cmd - t1');
sum5 = sumblk('u_Wt2perf = y_Wtdiffcmd + y_Wt1cmd - t2');

% This produces the uncertain state-space model
P = connect(tank1and2u,hot_act,cold_act,t1F,t2F,Wt1cmd,Wtdiffcmd,Whact, ...
                Wcact,Whrate,Wcrate,Wt1perf,Wt2perf,Wt1noise,Wt2noise, ...
                   sum1,sum2,sum3,sum4,sum5,inputs,outputs);

disp('Weighted open-loop model: ')
P

%% H-infinity Controller Design
nmeas = 4;		% Number of measurements
nctrls = 2;		% Number of controls
[k0,g0,gamma0] = hinfgsyn(P.NominalValue,nmeas,nctrls,0.901587304);
size(g0)
gamma0

inputs = {'t1ref', 't2ref', 't1noise', 't2noise', 'fhc', 'fcc'};
outputs = {'y_tank1', 'y_tank2', 'fhc', 'fcc', 'y_t1ref', 'y_t2ref', ...
                'y_t1Fn', 'y_t2Fn'};

hot_act(1).InputName = 'fhc'; hot_act(1).OutputName = 'y_hot_act';
cold_act(1).InputName = 'fcc'; cold_act(1).OutputName = 'y_cold_act';

tank1nom.InputName = [hot_act(1).OutputName cold_act(1).OutputName];
tank1nom.OutputName = 'y_tank1';
tank2nom.InputName = tank1nom.OutputName;
tank2nom.OutputName = 'y_tank2';

t1F.InputName = tank1nom.OutputName(2); t1F.OutputName = 'y_t1F';
t2F.InputName = tank2nom.OutputName; t2F.OutputName = 'y_t2F';

I_tref = zpk(eye(2));
I_tref.InputName = {'t1ref', 't2ref'}; I_tref.OutputName = {'y_t1ref', 'y_t2ref'};

sum1 = sumblk('y_t1Fn = y_t1F + t1noise');
sum2 = sumblk('y_t2Fn = y_t2F + t2noise');

simlft = connect(tank1nom,tank2nom,hot_act(1),cold_act(1),t1F,t2F,I_tref,sum1,sum2,inputs,outputs);

% Close the loop with the H-infinity controller |k0|
sim_k0 = lft(simlft,k0);
sim_k0.InputName = {'t1ref'; 't2ref'; 't1noise'; 't2noise'};
sim_k0.OutputName = {'h1'; 't1'; 't2'; 'fhc'; 'fcc'};

% Closed-loop simulation
time=0:800;
t1ref = (time>=80 & time<100).*(time-80)*-0.18/20 + ...
    (time>=100)*-0.18;
t2ref = (time>=80 & time<100).*(time-80)*-0.2/20 + ...
    (time>=100)*-0.2;
t1noise = Wt1noise.k * randn(size(time));
t2noise = Wt2noise.k * randn(size(time));

y = lsim(sim_k0,[t1ref ; t2ref ; t1noise ; t2noise],time);

h1 = h1ss+y(:,1);
t1 = t1ss+y(:,2);
t2 = t2ss+y(:,3);
fhc = fhss/fs+y(:,4); % Note scaling to actuator
fcc = fcss/fs+y(:,5); % Limits (0<= fhc <= 1) etc.

% figure
plot(time,h1,'--',time,t1,'-',time,t2,'-.');
xlabel('Time (sec)')
ylabel('Measurements')
title('Step Response of H-infinity Controller k0')
legend('h1','t1','t2');
grid

% figure
plot(time,fhc,'-',time,fcc,'-.');
xlabel('Time: seconds')
ylabel('Actuators')
title('Actuator Commands for H-infinity Controller k0')
legend('fhc','fcc');
grid

%% Robustness of the H-infinity Controller
clpk0 = lft(P,k0);

% Compute and plot worst-case gain
figure
wcsigma(clpk0);
axis([1e-4 100 -20 10]), shg