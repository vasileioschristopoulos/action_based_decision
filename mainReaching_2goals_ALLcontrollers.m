function mainReaching_2goals_ALLcontrollers(debug)

close all
clc

fig = figure('Position',[50,50,900,700],'Name','Neural Field',...
  'Color','w','NumberTitle','off','MenuBar','none');

% create axes for field plots
topAxes    = axes('Position',[0.1 0.65 0.85 0.3]);
bottomAxes = axes('Position',[0.1 0.3 0.85 0.3]);
trajAxes   = axes('Position',[0.1 0.3 0.85 0.3]);

topAxes    = axes('Position',[0.1 0.75 0.85 0.2]);
bottomAxes = axes('Position',[0.1 0.45  0.85 0.2]);

% create sliders for model parameters
controlFieldHeight = 0.04;
controlFieldWidth = 0.3;
sliderWidth = 0.2;
gapWidth = 0.01;
textWidth = controlFieldWidth - sliderWidth - gapWidth;

controlParamNames = {'dnf.params.h_u','dnf.params.c_exc','dnf.params.c_inh','dnf.params.g_inh','dnf.params.q_u','dnf.params.beta_u'};
controlPosX = [0, 0, 2, 2, 0, 0] * controlFieldWidth;
controlPosY = [5, 4, 5, 4, 3, 2] * controlFieldHeight;
controlMin = [-10, 0, 0, 0, 0, 0];
controlMax = [0, 100, 100, 5, 1.5, 5.0];
textFormat = {'%0.1f', '%0.1f', '%0.1f', '%0.2f', '%0.2f', '%0.2f'};

if debug==0
    close all
end

control_point = 0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% define parameters
N = 100;                % duration in number of time steps
prm.dt = 0.01;          % time step (sec)
prm.c  = 0.01;           % control-dependent noise
prm.r1  = 0.1;       % control signal penalty Ux
prm.r2  = 0.1;       % control signal penalty Uy

prm.Xh     = 0.5*0.02;  % position X noise
prm.Yh     = 0.5*0.02;  % position Y noise
prm.Vx     = 0.5*0.2;   % velocity X noise
prm.Vy     = 0.5*0.2;   % velocity Y noise
prm.frcX   = 0.5*1.0;   % force X noise
prm.frcY   = 0.5*1.0;   % force Y noise

NSim  = 1;            % Number of simulated trajectories
Init  = 1;            % LQG parametet (Initialization *1*)
Niter = 0;            % LQG parametet (Num of iterations *0*)

m     = 1;
tau_1 = 0.04;
tau_2 = 0.04;

%%Targets position and hand fixation
xh = 0;           % Hand initial x-poisition
yh = 0;           % Hand initial y-poisition

rx1 = -12;        % Initial x-position of the target 1
ry1 =  40;        % Initial y-position of the target 1

rx2 = 12;         % Initial x-position of the target 2
ry2 = 40;         % Initial y-position of the target 2

%Action-related decision variables
lambda   = 1;      % Inverse temperature variable (scaling the action cost)

%move_threshold      = 0.75;    % Threshold for start moving
move_threshold      = 0.6;    % Threshold for start moving
%move_threshold      = 0.95;    % Threshold for start moving
%move_threshold      = 90;    % Threshold for start moving
%output_u_threshold  = 0.5;   % Threshold for activating the associated controlles
output_u_threshold  = 0.25;   % Threshold for activating the associated controllers
%output_u_threshold  = 15;   % Threshold for activating the associatedcontrollers
revaluate_threshold = 10;     % Every n-time steps re-evaluate the policies
Threshold_contact   = 6;      % Threshold for contact with the target
flag_move           = 0.0;
cnt_evaluate        = 1.0;    % Counter for re-evaluation
flag_evaluate       = 1.0;
flag_contact        = 0.0;    % When 1 --> effector arrives to the target

%Parameters for the history of trials
learn_rate = 0.005;
Seq        = 5;



%%%%%%%%%%%%%%%%%%%%
% model parameters %
%%%%%%%%%%%%%%%%%%%%

dnfParams=initDNFParams();

%%%%%%%%%%%%%%%%%%%%%%%%%%
% simulation time course %
%%%%%%%%%%%%%%%%%%%%%%%%%%

nTrials = 1;
tMax = 300;

% set times at which field activities are stored (different variants)
% tStoreFields = [100, 200]; % select specific time steps
tStoreFields = 1:tMax; % store field activities at every time step
% tStoreFields = 1:5:tMax; % store field activities every 5th time step

% set times of stimulus presentation
% for multiple separate stimuli in one trial, repeat this for each one
% (e. g. tStimulusStart1 = ..., tStimulusStart2 = ...)
tStimulusStart = 50;
tStimulusEnd   = 400;
tStimulusDuration = tStimulusEnd-tStimulusStart;

%%%%%%%%%%%%%%%%%%
% initialization %
%%%%%%%%%%%%%%%%%%

dnf=initDNF(dnfParams, nTrials, tStoreFields);


%%%%%%%%%%%%%%%%%%%%
% Init Controllers %
%%%%%%%%%%%%%%%%%%%%
nControlParams = length(controlParamNames);
sliders = zeros(nControlParams, 1);
textFields = zeros(nControlParams, 1);
if debug>0
  for i = 1 : nControlParams
    eval(['tmp = ' controlParamNames{i} ';']);
    sliders(i) = uicontrol(fig, 'Style', 'Slider', 'Units', 'Norm', 'Position', ...
      [controlPosX(i)+textWidth+gapWidth, controlPosY(i), sliderWidth, controlFieldHeight], ...
      'Value', tmp, 'Min', controlMin(i), 'Max', controlMax(i), 'Callback', @sliderCallback);
    textFields(i) = uicontrol(fig,'Style','Text','Units','Norm','HorizontalAlignment', 'left', ...
      'String',[controlParamNames{i} '=' num2str(tmp, textFormat{i})], 'BackgroundColor', 'w',...
      'Position',[controlPosX(i)+gapWidth, controlPosY(i) textWidth controlFieldHeight]);
  end
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

%%%%%%%%%%%%%%
% simulation %
%%%%%%%%%%%%%%

Distance_from_origin_T1 = sqrt((X1(1)-rx1)^2  + (X1(2)-ry1)^2); %current distance from target 1
Distance_from_origin_T2 = sqrt((X1(1)-rx2)^2  + (X1(2)-ry2)^2); %current distance from target 1
distance_from_origin    = min(Distance_from_origin_T1,Distance_from_origin_T2);

% Protocol, 0=free choice, 1=instructed, 2=cued target uncertainty, 3=target uncertainty, 4=target prob, 5=target history
protocol=5;
stimulus=zeros(1,dnf.params.fieldSize);

% plot graphs
actPlot_u=0;
outPlot_u=0;
inPlot=0;
if debug>0
  axes(topAxes);
  cla;
  hold on;
  plot([0,dnf.params.fieldSize-1],[0,0],'Linestyle',':','Linewidth',1);
  actPlot_u = plot(0:dnf.params.fieldSize-1,dnf.field_u,'color','b','Linewidth',3);
  outPlot_u = plot(0:dnf.params.fieldSize-1,10*dnf.output_u,'color','r','Linewidth',1);
  inPlot = plot(0:dnf.params.fieldSize-1,stimulus+dnf.params.h_u,'color','g','Linewidth',1);
  set(gca,'ylim',[-15,15],'xlim',[0,dnf.params.fieldSize-1],'Ytick',[-10,0,10]);
  ylabel('u field','Fontsize',12);
  hold off;

  axes(bottomAxes);
  cla;
  hold on;
  plot([-dnf.params.halfField,dnf.params.halfField],[0,0],'Linestyle',':','Linewidth',1);
  kernelPlot = plot(-dnf.params.halfField:dnf.params.halfField, ...
      [zeros(1, dnf.params.halfField-dnf.kSize_uu) dnf.kernel_uu zeros(1, dnf.params.halfField-dnf.kSize_uu)] - dnf.params.g_inh, 'Color', 'r', 'Linewidth', 3);
  set(gca,'ylim',[-10,10],'xlim',[-dnf.params.halfField,dnf.params.halfField],'Ytick',[-10,-5,0,5,10]);
  ylabel('interaction kernel','Fontsize',12);
  hold off;
end

% loop over trials
for i = 1 : nTrials
    
    disp(['trial ' num2str(i)]);
    
    dnf=resetDNF(dnf);
    
    % State vector
    xh = 0;           % Hand initial x-poisition
    yh = 0;           % Hand initial y-poisition
    X1      = zeros(10,1);
    X1(1)   = xh;
    X1(2)   = yh;
    X1(9)   = rx1;
    X1(10)  = ry1;

    flag_move           = 0.0;
    cnt_evaluate        = 1.0;    % Counter for re-evaluation
    flag_evaluate       = 1.0;
    flag_contact        = 0.0;    % When 1 --> effector arrives to the target

    RelDesir   =zeros(1,dnf.params.fieldSize);
    stimulus  =0.0;

    stim1Coeff=5.0;
    stim2Coeff=5.0;
    cueStim1Coeff=12.0;
    cueStim2Coeff=12.0;
    biasStim1Coeff=5.0;
    biasStim2Coeff=5.0;
    if protocol==0
        move_threshold      = 0.8;
    elseif protocol==1
        move_threshold      = 0.8;
        if rand()<.5
            stim1Coeff=0.0;
            cueStim1Coeff=0.0;
        else
            stim2Coeff=0.0;
            cueStim2Coeff=0.0;
        end
    elseif protocol==2 || protocol==3
        if rand()<.5
            cueStim1Coeff=0.0;
        else
            cueStim2Coeff=0.0;
        end
    elseif protocol==4
        biasStim1Coeff=5;
        biasStim2Coeff=5.5;
        cueStim2Coeff=0.0;
    elseif protocol==5
        biasStim2Coeff = 5;
        biasStim1Coeff = UpdateTargetProb(biasStim1Coeff,biasStim2Coeff,learn_rate,Seq);
        cueStim2Coeff  =0.0;
    end
    
    clear Trajectory;

    % loop over time steps
    for t = 1 : tMax
        
        %Target related biases
        [theta1,rho1] = cart2pol(rx1-X1(1),ry1-X1(2));
        theta1=(theta1/pi)*180.0;
        [theta2,rho2] = cart2pol(rx2-X1(1),ry2-X1(2));
        theta2=(theta2/pi)*180.0;
        stim1 = stim1Coeff*gauss(1:dnf.params.fieldSize, round(theta1), 5); % a localized input
        stim2 = stim2Coeff*gauss(1:dnf.params.fieldSize, round(theta2), 5); % a localized input
        cue_stim1 = cueStim1Coeff*gauss(1:dnf.params.fieldSize, round(theta1), 5); % Cued stimulus bias for T1
        cue_stim2 = cueStim2Coeff*gauss(1:dnf.params.fieldSize, round(theta2), 5); % Cued stimulus bias for T2
        bias_stim1 = biasStim1Coeff*gauss(1:dnf.params.fieldSize, round(theta1), 5); % Biased of expected reward to T1
        bias_stim2 = biasStim2Coeff*gauss(1:dnf.params.fieldSize, round(theta2), 5); % Biased of expected reward to T2

        if t<tStimulusStart 
            stimulus   = bias_stim1+bias_stim2;
        elseif  t>= tStimulusStart && t<tStimulusEnd
            stimulus   = bias_stim1+bias_stim2+stim1+stim2;
            if (protocol==2 || protocol==4 || protocol==5) && (t-tStimulusStart)/tStimulusDuration>=.2
                stimulus = stimulus + cue_stim1+cue_stim2;
            end
            if protocol==3 && (t-tStimulusStart)/tStimulusDuration>=.3
                if cueStim1Coeff>0.0
                    stimulus=bias_stim1+bias_stim2+stim1+cue_stim1;
                else
                    stimulus=bias_stim1+bias_stim2+stim2+cue_stim2;
                end
            end
            RelDesirField = 0;
            for myFieldSize = 1:dnf.params.fieldSize
                theta     = myFieldSize;
                RelDesirField = RelDesirField + .25*RelDesir(myFieldSize)*gauss(1:dnf.params.fieldSize, round(theta), 5); % a localized input
            end
            stimulus     = stimulus + RelDesirField;
        elseif t>tStimulusEnd
            stimulus = zeros(1,dnf.params.fieldSize);
        end
        
        
        disp(['t=' num2str(t)]);
        dnf=runDNF(dnf, stimulus, tStoreFields, t);
    
        if debug>0
            set(inPlot, 'Ydata', stimulus+dnf.params.h_u);
            set(actPlot_u,'Ydata',dnf.field_u);
            set(outPlot_u,'Ydata',10*dnf.output_u);
            drawnow;
        end
        %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %LQG controller
        if t>=tStimulusStart && flag_contact == 0
            Cost_commands = 0; Cost_accuracy =0; Overall_action_cost =0;
            if max(dnf.output_u)<move_threshold  & flag_move ==0   %effector has not moved yet (planning of actions)
                u_above_thr_planning_id = find(dnf.output_u>output_u_threshold);
                if isempty(u_above_thr_planning_id) == 0 %Planning
                    for myPlanning = 1:length(u_above_thr_planning_id)
                        [x y] = pol2cart(u_above_thr_planning_id(myPlanning)*pi/180, distance_from_origin);
                        x=x+X1(1); y=y+X1(2);
                        X1(end-1:end) = [x y];
                        [K,L,Cost,Xa,XSim,Xhat,CostSim,Unoise,Ufree] = kalman_lqg(A,B,C,C0, H,D,D0, E0, Q,R, X1,S1, NSim,Init,Niter );
                        [Cost_commands(myPlanning) Cost_accuracy(myPlanning) Overall_action_cost(myPlanning)] = ActionCost(Unoise,Q,R,Xhat,[rx1 ry1],[rx2 ry2]);
                    end
                    Overall_action_cost                = exp(-Overall_action_cost);
                    NormalFactor                       = sum(Overall_action_cost);
                    RelDesir                           = zeros(1,dnf.params.fieldSize);
                    RelDesir(u_above_thr_planning_id)  = Overall_action_cost./NormalFactor;
                end
                
                %elseif max(dnf.output_u)>move_threshold %effector starts moving
            else %effector starts moving
                %move
                flag_move = 1;
                if flag_evaluate          == 1 % Re-evaluate the policies
                    cnt_evaluate           = 1;
                    flag_evaluate          = 0;
                    u_above_thr_id         = find(dnf.output_u>output_u_threshold);
                    u_above_thr            = dnf.output_u(u_above_thr_id)/sum(dnf.output_u(u_above_thr_id));
                    NumOfActiveControllers = length(u_above_thr_id);
                    U_avg                  = 0;  %weighted average of the motor commands
                    current_distance_T1    = sqrt((X1(1)-rx1)^2  + (X1(2)-ry1)^2); %current distance from target 1
                    current_distance_T2    = sqrt((X1(1)-rx2)^2  + (X1(2)-ry2)^2); %current distance from target 2
                    Run_distance_t1(t)        = current_distance_T1;%min(current_distance_T1,current_distance_T2);
                    Run_distance_t2(t)        = current_distance_T2;
                    disp(['Distance_to_contact=' num2str(min(Run_distance_t1(t),Run_distance_t2(t)))]);
                    if Run_distance_t1(t) < Threshold_contact || Run_distance_t2(t) < Threshold_contact
                        flag_contact = 1;
                        Contact_time = t;
                        tStimulusEnd = Contact_time; %When we contact the target, the stimuli disappear
                    end
                    for myControl = 1:NumOfActiveControllers
                        [x y] = pol2cart(u_above_thr_id(myControl)*pi/180, max(Run_distance_t1(t),Run_distance_t2(t)));
                        x=x+X1(1); y=y+X1(2);
                        X1(end-1:end) = [x y];
                        [K,L,Cost,Xa,XSim,Xhat,CostSim,Unoise,Ufree] = kalman_lqg(A,B,C,C0, H,D,D0, E0, Q,R, X1,S1, NSim,Init,Niter );
                        [Cost_commands(myControl) Cost_accuracy(myControl) Overall_action_cost(myControl)] = ActionCost(Unoise,Q,R,Xhat,[rx1 ry1],[rx2 ry2]);
                        U_avg   = U_avg +  u_above_thr(myControl)*squeeze(Unoise);
                    end
                    if NumOfActiveControllers>0
                        Xgo                                = SimTraj(A,B,C0,H,Q,U_avg,NSim,X1);
                        Trajectory(:,t)                    = Xgo(1:2,cnt_evaluate);
                        Overall_action_cost                = exp(-Overall_action_cost);
                        NormalFactor                       = sum(Overall_action_cost);
                        RelDesir                           = zeros(1,dnf.params.fieldSize);
                        RelDesir(u_above_thr_id)           = Overall_action_cost./NormalFactor;
                    else
                        Trajectory(:,t)                    = Trajectory(:,t-1);
                    end
                else
                    cnt_evaluate          = cnt_evaluate + 1;
                    Trajectory(:,t)       = Xgo(1:2,cnt_evaluate);
                    X1(1:8)               = Xgo(:,cnt_evaluate);     %Update the state vector of T1
                    current_distance_T1   = sqrt((X1(1)-rx1)^2  + (X1(2)-ry1)^2); %current distance from target 1
                    current_distance_T2   = sqrt((X1(1)-rx2)^2  + (X1(2)-ry2)^2); %current distance from target 2
                    Run_distance_t1(t)       = current_distance_T1;%min(current_distance_T1,current_distance_T2);
                    Run_distance_t2(t)       = current_distance_T2;
                    if Run_distance_t1(t) < Threshold_contact | Run_distance_t2(t) < Threshold_contact | (Run_distance_t1(t) > Run_distance_t1(t-1) && Run_distance_t2(t) > Run_distance_t2(t-1))
                        flag_contact = 1;
                        Contact_time = t;
                        tStimulusEnd = Contact_time;   %When we contact the target, the stimuli disappear
                        disp(['Distance_to_contact=' num2str(min(Run_distance_t1(t),Run_distance_t2(t)))]);
                    end
                    if cnt_evaluate      == revaluate_threshold      %Re-evaluate the policies on the next step
                        flag_evaluate      = 1;
                    end
                end
            end
        end
        
    end
    
    %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    file_name=['/home/jbonaiuto/Projects/actionBasedDecision/trial.' num2str(i) '.trajectory.mat'];
    save(file_name, 'Trajectory');
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
imagesc(dnf.history_output');
xlabel('time');
ylabel('output u');
colorbar();

% View plots of all stored field activities iteratively
if 0
    
end

% view evolution of field activities in each trial as mesh plot
nFieldsPerTrial = length(tStoreFields);
if 0
    disp('Press any key to iterate through trials');
    figure;
    for i = 1 : nStoredFields
        plot(0:dnf.params.fieldSize-1, zeros(1, dnf.params.fieldSize), ':k', ...
            1:dnf.params.fieldSize, dnf.history_s(i, :), '--g', ...
            1:dnf.params.fieldSize, dnf.history_mu(i, :) + dnf.params.h_u, '-c', ...
            1:dnf.params.fieldSize, dnf.history_u(i, :), '-b');
        set(gca, 'XLim', [0 dnf.params.fieldSize-1], 'YLim', [-15 15]);
        ylabel('activity u');
        drawnow;
        pause(0.01);
    end
    for i = 1 : nTrials
        subplot(2, 1, 1);
        mesh(1:dnf.params.fieldSize, tStoreFields, ...
            dnf.history_u((i-1)*nFieldsPerTrial+1 : i*nFieldsPerTrial, :));
        zlabel('activity u');
        subplot(2, 1, 2);
        mesh(1:dnf.params.fieldSize, tStoreFields, ...
            dnf.history_mu((i-1)*nFieldsPerTrial+1 : i*nFieldsPerTrial, :));
        zlabel('memory trace u');
        pause
    end
end

% create  movie of the activation u and trajectory
if 0
    figure;
    for i = 1 : nStoredFields
        h1 = subplot(2,1,1);
        plot(0:dnf.params.fieldSize-1, zeros(1, dnf.params.fieldSize), ':k', ...
            1:dnf.params.fieldSize, dnf.history_s(i, :), '--g', ...
            1:dnf.params.fieldSize, dnf.history_mu(i, :) + dnf.params.h_u, '-c', ...
            1:dnf.params.fieldSize, dnf.history_u(i, :), '-b');
        set(gca, 'XLim', [0 dnf.params.fieldSize-1], 'YLim', [-15 15]);
        ylabel('activity u');
        drawnow;
        pause(0.01);
        h2 = subplot(2,1,2);
        if i>= tStimulusStart && i<=Contact_time
            plot(0,0,'sr'); hold on
            circle([rx1 ry1],5,1000,3)
            circle([rx2 ry2],5,1000,3)
            plot(Trajectory(1,i),Trajectory(2,i),'ko','MarkerSize',3)
            axis equal
            pause(0.01);
        elseif i>Contact_time
            cla(h2)
            h2 = subplot(2,1,2);
            plot(0,0,'sr'); hold on
            plot(rx1,ry1,'ko','MarkerSize',3);
            plot(rx2,ry2,'ko','MarkerSize',3);
        end
    end
end

% view mesh plot of all stored field activities together
if 0
    figure;
    subplot(2, 1, 1);
    mesh(1:dnf.params.fieldSize, 1:nStoredFields, dnf.history_u(:, :));
    zlabel('activity u');
    subplot(2, 1, 2);
    mesh(1:dnf.params.fieldSize, 1:nStoredFields, dnf.history_mu(:, :));
    zlabel('memory trace u');
end

figure; hold on
circle([rx1 ry1],5,1000,3)
circle([rx2 ry2],5,1000,3)
plot(0,0,'sk')
plot(Trajectory(1,:),Trajectory(2,:),'k','LineWidth',3)
axis equal


% update paramter values after slider changed
function sliderCallback(hObject, eventdata) %#ok<INUSD>
  sliderChanged = find(hObject == sliders);

  paramName = controlParamNames{sliderChanged};
  tmp = get(sliders(sliderChanged), 'Value');
  set(textFields(sliderChanged), 'String', [paramName '=' num2str(tmp, textFormat{sliderChanged})]);
  eval([paramName '= tmp;']);
  dnf.kernel_uu = dnf.params.c_exc * gaussNorm(-dnf.params.halfField:dnf.params.halfField, 0, dnf.params.sigma_exc) ...
    - dnf.params.c_inh * gaussNorm(-dnf.params.halfField:dnf.params.halfField, 0, dnf.params.sigma_inh) - dnf.params.g_inh;
  dnf.kernel_mu = dnf.params.c_mu * gaussNorm(-dnf.params.halfField:dnf.params.halfField, 0, dnf.params.sigma_mu);
  set(kernelPlot,'YData',[zeros(1, dnf.params.halfField-dnf.kSize_uu) dnf.kernel_uu zeros(1, dnf.params.halfField-dnf.kSize_uu)] - dnf.params.g_inh);
end
end



