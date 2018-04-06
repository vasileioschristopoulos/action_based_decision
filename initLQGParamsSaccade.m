function lqgParams_saccade=initLQGParamsSaccade()

N = 100;               % duration in number of time steps
prm.dt = 0.001;        % time step (sec)
prm.c  = 0.01;          % control-dependent noise
prm.r1  = 0.001;         % control signal penalty torque x
prm.r2  = 0.001;         % control signal penalty torque y

prm.Xe     = 0.00001;      % control signal penalty Xe
prm.Ye     = 0.00001;      % control signal penalty Ye
prm.Vxe    = 0.00002;      % velocity signal penalty Xe vel
prm.Vye    = 0.00002;      % velocity signal penalty Ye vel

w_ac   = 1;
w_stop = 1;


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Eye plant model

%%----HUMAN-----
tau1 = 0.224;
tau2 = 0.013;
%%--------------

% %%---MONKEY-----
% tau1 = 0.260;
% tau2 = 0.012;
% %%--------------

kc = 1;
bc = tau1 + tau2;
mc = tau1*tau2;

lqgParams_saccade.NSim  = 1;   % Number of simulated trajectories
lqgParams_saccade.Init  = 1;   % LQG parametet (Initialization *1*)
lqgParams_saccade.Niter = 0;   % LQG parametet (Num of iterations *0*)


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Compute system dynamics and cost matrices

%Ac: Continuous, A: Discrete 
Ac       = zeros(4,4);
Ac(1,2)  = 1;
Ac(2,1)  = -kc/mc;
Ac(2,2)  = -bc/mc;
Ac(3,4)  = 1;
Ac(4,3)  = -kc/mc;
Ac(4,4)  = -bc/mc;

lqgParams_saccade.A = expm(Ac*prm.dt);
lqgParams_saccade.A = [lqgParams_saccade.A zeros(4,2);zeros(2,4) eye(2)];

%Bc: Continuous, B: Discrete
Bc       = zeros(4,2);
Bc(2,1)  = 1/mc;
Bc(4,2)  = 1/mc;

lqgParams_saccade.B  = inv(Ac)*(expm(Ac*prm.dt) - eye(4))*Bc;
lqgParams_saccade.B  = [lqgParams_saccade.B; zeros(2,2)];

lqgParams_saccade.C   = prm.c;
lqgParams_saccade.C0  = 0;

lqgParams_saccade.D = 0;
lqgParams_saccade.D0 = diag([prm.Xe prm.Ye prm.Vxe prm.Vye]);


lqgParams_saccade.H = zeros(4,6);
lqgParams_saccade.H(1:4,1:4) = eye(4);


lqgParams_saccade.E0 = 0;

lqgParams_saccade.R = diag([prm.r1,prm.r2])/N;


lqgParams_saccade.Q = zeros(6,6,N);

lqgParams_saccade.Q(:,:,N) = eye(6);
lqgParams_saccade.Q(1,1,N) = w_ac;
lqgParams_saccade.Q(2,2,N) = w_stop;
lqgParams_saccade.Q(3,3,N) = w_ac;
lqgParams_saccade.Q(4,4,N) = w_stop;
lqgParams_saccade.Q(5,5,N) = w_ac;
lqgParams_saccade.Q(6,6,N) = w_ac;


lqgParams_saccade.Q(1,5,N) = -w_ac;
lqgParams_saccade.Q(5,1,N) = -w_ac;
lqgParams_saccade.Q(3,6,N) = -w_ac;
lqgParams_saccade.Q(6,3,N) = -w_ac;


lqgParams_saccade.S1 = zeros(6,6);

