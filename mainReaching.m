clear all
close all
%clc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% define parameters
N = 250;                % duration in number of time steps
prm.dt = 0.01;          % time step (sec)
prm.c  = 1.0;           % control-dependent noise
prm.r1  = 0.0001;       % control signal penalty Ux
prm.r2  = 0.0001;       % control signal penalty Uy

prm.Xh     = 0.5*0.02;  % position X noise
prm.Yh     = 0.5*0.02;  % position Y noise
prm.Vx     = 0.5*0.2;   % velocity X noise
prm.Vy     = 0.5*0.2;   % velocity Y noise
prm.frcX   = 0.5*1.0;   % force X noise
prm.frcY   = 0.5*1.0;   % force Y noise

%NSim  = 100;          % Number of simulated trajectories
NSim  = 10;          % Number of simulated trajectories
Init  = 1;            % LQG parametet (Initialization *1*)
Niter = 0;            % LQG parametet (Num of iterations *0*)

m     = 1;
tau_1 = 0.04;
tau_2 = 0.04;

%%Targets position and hand fixation
xh = 0;           % Hand initial x-poisition
yh = 0;           % Hand initial y-poisition

rx1 = -20;         % Initial x-position of the target 1
ry1 =  40;        % Initial y-position of the target 1

rx2 = 20;         % Initial x-position of the target 2
ry2 =  40;        % Initial y-position of the target 2

%%%%%%%%%%%%%%%%%%%%
% model parameters %
%%%%%%%%%%%%%%%%%%%%

neuralParam.fieldSize = 181; % must be odd

neuralParam.tau_u = 20; % time constant of dynamic field
neuralParam.tau_mu_build = 500; neuralParam.tau_mu_decay = 2000; % time constants of memory traces
neuralParam.beta_u = 5; % steepness parameter of sigmoid function
neuralParam.h_u = -5; % resting level

% interation paramters
% c: strength of interaction kernel
% sigma: width of interaction kernel
% g: strength of global inhibition
neuralParam.c_exc = 15; neuralParam.sigma_exc = 5;
neuralParam.c_inh = 0; neuralParam.sigma_inh = 10; neuralParam.g_inh = 0.5;
neuralParam.c_mu = 1; neuralParam.sigma_mu = 5;

neuralParam.q_u = 0; % noise levels
neuralParam.sigma_q = 0; % width of the noise kernel

%%%%%%%%%%%%%%%%%%%%%%%%%%
% simulation time course %
%%%%%%%%%%%%%%%%%%%%%%%%%%

nTrials = 1;
tMax = 250;
tMax = 5;

% set times at which field activities are stored (different variants)
% tStoreFields = [100, 200]; % select specific time steps
tStoreFields = 1:tMax; % store field activities at every time step
% tStoreFields = 1:5:tMax; % store field activities every 5th time step

% set times of stimulus presentation
% for multiple separate stimuli in one trial, repeat this for each one
% (e. g. tStimulusStart1 = ..., tStimulusStart2 = ...)
tStimulusStart = 50;
tStimulusEnd = 150;

[theta1,rho1] = cart2pol(rx1-xh,ry1-yh);
theta1=(theta1/pi)*180.0;
stim1 = 6*gauss(1:neuralParam.fieldSize, round(theta1), 5); % a localized input
[theta2,rho2] = cart2pol(rx2-xh,ry2-yh);
theta2=(theta2/pi)*180.0;
stim2 = 6*gauss(1:neuralParam.fieldSize, round(theta2), 5); % a localized input

%%%%%%%%%%%%%%%%%%
% initialization %
%%%%%%%%%%%%%%%%%%

neuralParam.halfField = floor(neuralParam.fieldSize/2);

% create row vectors for field activities
dnf.field_u = zeros(1, neuralParam.fieldSize);
dnf.memTrace_u = zeros(1, neuralParam.fieldSize);

% create matrices to store field activities at different times
dnf.history_s = zeros(nTrials * length(tStoreFields), neuralParam.fieldSize);
dnf.history_u = zeros(nTrials * length(tStoreFields), neuralParam.fieldSize);
dnf.history_mu = zeros(nTrials * length(tStoreFields), neuralParam.fieldSize);

% index of the current position in the history matrices
dnf.iHistory = 1;

% set up the interaction kernels
dnf.kernel_uu = neuralParam.c_exc * gaussNorm(-neuralParam.halfField:neuralParam.halfField, 0, neuralParam.sigma_exc) ...
  - neuralParam.c_inh * gaussNorm(-neuralParam.halfField:neuralParam.halfField, 0, neuralParam.sigma_inh) - neuralParam.g_inh;
dnf.kernel_mu = neuralParam.c_mu * gaussNorm(-neuralParam.halfField:neuralParam.halfField, 0, neuralParam.sigma_mu);

% set up the kernel for correlated noise (if required)
if neuralParam.sigma_q > 0
  dnf.kernel_q = gaussNorm(-neuralParam.halfField:neuralParam.halfField, 0, neuralParam.sigma_q);
end

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% compute system dynamics and cost matrices
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

A = zeros(10,10);
A(1,1) = 1; A(1,3) = prm.dt;
A(2,2) = 1; A(2,4) = prm.dt;
A(3,3) = 1; A(3,5) = prm.dt/m;
A(4,4) = 1; A(4,6) = prm.dt/m;
A(5,5) = 1-(prm.dt/tau_2);  A(5,7) = prm.dt/tau_2;
A(6,6) = 1-(prm.dt/tau_2);  A(6,8) = prm.dt/tau_2;
A(7,7) = 1-(prm.dt/tau_1);
A(8,8) = 1-(prm.dt/tau_1);
A(9,9) = 1;
A(10,10) = 1;


B = zeros(10,2);
B(7,1) = prm.dt/tau_1;
B(8,2) = prm.dt/tau_2;

C = prm.c;
C0 = 0;


H = zeros(6,10);
H(1:6,1:6) = eye(6);


D  = 0;
D0 = diag([prm.Xh prm.Yh prm.Vx prm.Vy prm.frcX prm.frcY]);

E0 = 0;

R = diag([prm.r1,prm.r2])/N;


Q = zeros(10,10,N);

Q(:,:,N) = eye(10);


Q(1,9,N) = -1;
Q(2,10,N) = -1;
Q(9,1,N) = -1;
Q(10,2,N) = -1;
Q(3,3,N) = 1;
Q(4,4,N) = 1;


S1 = zeros(10,10);

% State vector 
X1      = zeros(10,1);
X1(1)   = xh;
X1(2)   = yh;
X1(9)   = rx1;
X1(10)  = ry1;

X2      = zeros(10,1);
X2(1)   = xh;
X2(2)   = yh;
X2(9)   = rx1;
X2(10)  = ry1;

%%%%%%%%%%%%%%
% simulation %
%%%%%%%%%%%%%%

% loop over trials
for i = 1 : nTrials
  
  disp(['trial ' num2str(i)]);

  % prepare matrix that holds the stimulus for each time step
  stimulus = zeros(tMax, neuralParam.fieldSize);
  
  % write the stimulus pattern into the stimulus matrix for all time steps
  % where it should be active
  for j = tStimulusStart : tStimulusEnd
    stimulus(j, :) = stim1+stim2;
  end
  
  % reset field activities to resting levels
  dnf.field_u(1:neuralParam.fieldSize) = neuralParam.h_u;
  
  % loop over time steps
  for t = 1 : tMax
    disp(['t=' num2str(t)]);
    % calculation of field outputs
    dnf.output_u = sigmoid(dnf.field_u, neuralParam.beta_u, 0);
    
    % circular padding of outputs for convolution
    dnf.output_u_padded = [dnf.output_u(neuralParam.halfField+2:neuralParam.fieldSize), dnf.output_u, ...
        dnf.output_u(:, 1:neuralParam.halfField)];
    dnf.memTrace_u_padded = [dnf.memTrace_u(neuralParam.halfField+2:neuralParam.fieldSize), dnf.memTrace_u, ...
        dnf.memTrace_u(:, 1:neuralParam.halfField)];

    % get endogenous input to fields by convolving outputs with interaction kernels
    dnf.conv_uu = conv2(1, dnf.kernel_uu, dnf.output_u_padded, 'valid');
    dnf.conv_mu = conv2(1, dnf.kernel_mu, dnf.memTrace_u_padded, 'valid');
        
    % create field noise for this timestep
    dnf.noise_u = neuralParam.q_u * randn(1, neuralParam.fieldSize);
    if neuralParam.sigma_q > 0 % create spatially correlated noise by convolution
      dnf.noise_u_padded = [dnf.noise_u(neuralParam.halfField+2:neuralParam.fieldSize), dnf.noise_u, ...
         dnf.noise_u(:, 1:neuralParam.halfField)];
      dnf.noise_u = conv2(1, dnf.kernel_q, dnf.noise_u_padded, 'valid');
    end
    
    % update field activities
    dnf.field_u = dnf.field_u + 1/neuralParam.tau_u * (-dnf.field_u + neuralParam.h_u + stimulus(t, :) + ...
        dnf.conv_uu + dnf.conv_mu) + dnf.noise_u;
        
    % update memory trace (only if there is activity in the field)
    dnf.activeRegions_u = dnf.field_u > 0;
    if any(dnf.activeRegions_u)
      dnf.memTrace_u = dnf.memTrace_u + 1/neuralParam.tau_mu_build * (-dnf.memTrace_u + dnf.output_u) .* ...
          dnf.activeRegions_u + 1/neuralParam.tau_mu_decay * (-dnf.memTrace_u) .* (1-dnf.activeRegions_u);
    end
    
    % store field activities at the selected time steps
    if any(tStoreFields == t)
      dnf.history_s(dnf.iHistory, :) = stimulus(t, :);
      dnf.history_u(dnf.iHistory, :) = dnf.field_u;
      dnf.history_mu(dnf.iHistory, :) = dnf.memTrace_u;
      dnf.iHistory = dnf.iHistory + 1;
    end

    %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %LQG controller
    
    %optimal control target 1 (T1)
    [K_T1,L_T1,Cost_T1,Xa_T1,XSim_T1,Xhat_T1,CostSim_T1,Unoise_T1,Ufree_T1] = ...
        kalman_lqg(A,B,C,C0, H,D,D0, E0, Q,R, X1,S1, NSim,Init,Niter );

    [Cost_commands_T1 Cost_accuracy_T1] = ActionCost(Unoise_T1,Q,R,X1);
    OverallCost_T1 = Cost_commands_T1 + Cost_accuracy_T1;
    
    %optimal control target 2 (T2)
    [K_T2,L_T2,Cost_T2,Xa_T2,XSim_T2,Xhat_T2,CostSim_T2,Unoise_T2,Ufree_T2] = ...
        kalman_lqg(A,B,C,C0, H,D,D0, E0, Q,R, X2,S1, NSim,Init,Niter );
    
    [Cost_commands_T2 Cost_accuracy_T2] = ActionCost(Unoise_T2,Q,R,X2);
    OverallCost_T2 = Cost_commands_T2 + Cost_accuracy_T2;
    
    RelActionCostT1 = exp(-OverallCost_T1/sigma)/(exp(-OverallCost_T1/sigma) + exp(-OverallCost_T2/sigma));
    RelActionCostT2 = exp(-OverallCost_T2/sigma)/(exp(-OverallCost_T1/sigma) + exp(-OverallCost_T2/sigma));

    %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  end
end
nStoredFields = nTrials * length(tStoreFields);

figure;
subplot(3, 1, 1);
imagesc(dnf.history_s');
xlabel('time');
ylabel('stimulus');
colorbar();
subplot(3, 1, 2);
imagesc(dnf.history_u');
xlabel('time');
ylabel('activity u');
colorbar();
subplot(3, 1, 3);
imagesc(dnf.history_mu');
xlabel('time');
ylabel('mem trace u');
colorbar();

% View plots of all stored field activities iteratively
if 0
  figure;
  for i = 1 : nStoredFields
    plot(0:neuralParam.fieldSize-1, zeros(1, neuralParam.fieldSize), ':k', ...
      1:neuralParam.fieldSize, dnf.history_s(i, :), '--g', ...
      1:neuralParam.fieldSize, dnf.history_mu(i, :) + neuralParam.h_u, '-c', ...
      1:neuralParam.fieldSize, dnf.history_u(i, :), '-b');
    set(gca, 'XLim', [0 neuralParam.fieldSize-1], 'YLim', [-15 15]);
    ylabel('activity u');
    drawnow;
    pause(0.01);
  end
end

% view evolution of field activities in each trial as mesh plot
nFieldsPerTrial = length(tStoreFields);
if 0
  disp('Press any key to iterate through trials');
  figure;
  for i = 1 : nTrials
    subplot(2, 1, 1);
    mesh(1:neuralParam.fieldSize, tStoreFields, ...
        dnf.history_u((i-1)*nFieldsPerTrial+1 : i*nFieldsPerTrial, :));
    zlabel('activity u');
    subplot(2, 1, 2);
    mesh(1:neuralParam.fieldSize, tStoreFields, ...
        dnf.history_mu((i-1)*nFieldsPerTrial+1 : i*nFieldsPerTrial, :));
    zlabel('memory trace u');
    pause
  end
end

% view mesh plot of all stored field activities together
if 0
  figure;
  subplot(2, 1, 1);
  mesh(1:neuralParam.fieldSize, 1:nStoredFields, dnf.history_u(:, :));
  zlabel('activity u');
  subplot(2, 1, 2);
  mesh(1:neuralParam.fieldSize, 1:nStoredFields, dnf.history_mu(:, :));
  zlabel('memory trace u');
end

%figure
%for mySimTraj = 1: NSim
%    plot(squeeze(Xhat(1,mySimTraj,:)),squeeze(Xhat(2,mySimTraj,:))); hold on  %Simulated trajectories
%end
%plot(Xa(1,:),Xa(2,:),'r','LineWidth',3) % Average trajectory





