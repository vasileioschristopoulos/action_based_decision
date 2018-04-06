function lqgParams_reach=initLQGParamsReach()

N = 100;                % duration in number of time steps
prm.dt = 0.005;          % time step (sec)
prm.c  = 1.5;          % control-dependent noise
prm.r1  = 0.001;          % control signal penalty Ux
prm.r2  = 0.001;          % control signal penalty Uy

prm.Xh     = 0.5*0.02;  % position X noise
prm.Yh     = 0.5*0.02;  % position Y noise
prm.Vx     = 0.5*0.2;   % velocity X noise
prm.Vy     = 0.5*0.2;   % velocity Y noise
prm.frcX   = 0.5*1.0;   % force X noise
prm.frcY   = 0.5*1.0;   % force Y noise

lqgParams_reach.NSim  = 1;   % Number of simulated trajectories
lqgParams_reach.Init  = 1;   % LQG parametet (Initialization *1*)
lqgParams_reach.Niter = 0;   % LQG parametet (Num of iterations *0*)

m     = 1;
tau_1 = 0.04;
tau_2 = 0.04;

lqgParams_reach.A = zeros(10,10); 
lqgParams_reach.A(1,1) = 1; 
lqgParams_reach.A(1,3) = prm.dt; 
lqgParams_reach.A(2,2) = 1; 
lqgParams_reach.A(2,4) = prm.dt;
lqgParams_reach.A(3,3) = 1; 
lqgParams_reach.A(3,5) = prm.dt/m; 
lqgParams_reach.A(4,4) = 1; 
lqgParams_reach.A(4,6) = prm.dt/m; 
lqgParams_reach.A(5,5) = 1-(prm.dt/tau_2);  
lqgParams_reach.A(5,7) = prm.dt/tau_2; 
lqgParams_reach.A(6,6) = 1-(prm.dt/tau_2);
lqgParams_reach.A(6,8) = prm.dt/tau_2; 
lqgParams_reach.A(7,7) = 1-(prm.dt/tau_1); 
lqgParams_reach.A(8,8) = 1-(prm.dt/tau_1);
lqgParams_reach.A(9,9) = 1; 
lqgParams_reach.A(10,10) = 1;


lqgParams_reach.B = zeros(10,2); 
lqgParams_reach.B(7,1) = prm.dt/tau_1; 
lqgParams_reach.B(8,2) = prm.dt/tau_2;
lqgParams_reach.C = prm.c; 
lqgParams_reach.C0 = 0; 
lqgParams_reach.H = zeros(6,10); 
lqgParams_reach.H(1:6,1:6) = eye(6);
lqgParams_reach.D  = 0; 
lqgParams_reach.D0 = diag([prm.Xh prm.Yh prm.Vx prm.Vy prm.frcX prm.frcY]);
lqgParams_reach.E0 = 0; 
lqgParams_reach.R = diag([prm.r1,prm.r2])/N;

lqgParams_reach.Q = zeros(10,10,N); 
lqgParams_reach.Q(:,:,N) = eye(10); 
lqgParams_reach.Q(1,9,N) = -1; 
lqgParams_reach.Q(2,10,N) = -1;
lqgParams_reach.Q(9,1,N) = -1; 
lqgParams_reach.Q(10,2,N) = -1; 
lqgParams_reach.Q(3,3,N) = 1; 
lqgParams_reach.Q(4,4,N) = 1;

lqgParams_reach.S1 = zeros(10,10);
