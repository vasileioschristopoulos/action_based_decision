function MultipleDNFs_reward(debug, init_cue_w)

close all
clc

param_fig = figure('Position',[50,10,1400,300],'Name','DNF Params',...
    'Color','w','NumberTitle','off','MenuBar','none');

fig = figure('Position',[50,50,1200,700],'Name','Neural Field',...
    'Color','w','NumberTitle','off','MenuBar','none');

% create axes for field plots
SpatialInAxes_Sac = axes('Position',[0.05 0.88 0.19 0.1]);
SpatialInAxes_Rch = axes('Position',[0.27 0.88 0.18 0.1]);
SacAxes       = axes('Position',[0.05 0.73 0.19 0.1]);
ReachAxes     = axes('Position',[0.27 0.73 0.18 0.1]);
EffortSacAxes = axes('Position',[0.05 0.58 0.19 0.1]);
EffortRchAxes = axes('Position',[0.27 0.58 0.19 0.1]);
RewardAlloAxes= axes('Position',[0.05 0.43 0.09 0.1]);
RewardBodyAxes = axes('Position',[0.16 0.43 0.09 0.1]);
RewardEgoAxes = axes('Position',[0.27 0.43 0.19 0.1]);
kernelAxes    = axes('Position',[0.05 0.28 0.40 0.1]);
CueAxes    = axes('Position',[0.50 0.88 0.20 0.1]);
CueWAxes      = axes('Position',[0.75 0.88 0.20 0.1]);
BehaviorAxes  = axes('Position',[0.50 0.28 0.45 0.5]);


% %create sliders for model parameters
controlFieldHeight = 0.04;
controlFieldWidth = 0.25;
sliderHeight = .1;
sliderWidth = 0.1;
gapWidth = 0.01;
textWidth = controlFieldWidth - sliderWidth - gapWidth;

controlParamNames = {'dnf_in_e.params.h_u','dnf_in_e.params.c_exc','dnf_in_e.params.c_inh','dnf_in_e.params.g_inh','dnf_in_e.params.q_u',...    
    'dnf_in_e.params.beta_u','dnf_in_h.params.h_u','dnf_in_h.params.c_exc','dnf_in_h.params.c_inh','dnf_in_h.params.g_inh',...
    'dnf_in_h.params.q_u','dnf_in_h.params.beta_u','dnf_sac.params.h_u','dnf_sac.params.c_exc','dnf_sac.params.c_inh','dnf_sac.params.g_inh',...
    'dnf_sac.params.q_u','dnf_sac.params.beta_u','dnf_rch.params.h_u','dnf_rch.params.c_exc','dnf_rch.params.c_inh','dnf_rch.params.g_inh','dnf_rch.params.q_u',...
    'dnf_rch.params.beta_u'};
controlPosX = [0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 3, 3, 3, 3, 3, 3] * controlFieldWidth;
controlPosY = [8, 7, 6, 5, 4, 3, 8, 7, 6, 5, 4, 3, 8, 7, 6, 5, 4, 3, 8, 7, 6, 5, 4, 3] * sliderHeight;
controlMin = [-10, 0, 0, 0, 0, 0, -10, 0, 0, 0, 0, 0, -10, 0, 0, 0, 0, 0, -10, 0, 0, 0, 0, 0];
controlMax = [0, 100, 100, 5, 1.5, 5.0, 0, 100, 100, 5, 1.5, 5.0, 0, 100, 100, 5, 1.5, 5.0, 0, 100, 100, 5, 1.5, 5.0];
textFormat = {'%0.1f', '%0.1f', '%0.1f', '%0.2f', '%0.2f', '%0.2f','%0.1f', '%0.1f', '%0.1f', '%0.2f', '%0.2f', '%0.2f','%0.1f', '%0.1f', '%0.1f', '%0.2f', '%0.2f', '%0.2f','%0.1f', '%0.1f', '%0.1f', '%0.2f', '%0.2f', '%0.2f'};

if debug==0
    close all
end

%%%%%%%%%%%%%%%%%%%%
% model parameters %
%%%%%%%%%%%%%%%%%%%%

dnf_inputParams=initDNFParams();  %Spatial input field
dnf_inputParams.tau_u=5;
dnf_inputParams.c_inh = 0;
dnf_inputParams.c_exc = 0;

dnf_expRewParams=initDNFParams();
dnf_expRewParams.tau_u=5;
dnf_expRewParams.c_inh = 0;
dnf_expRewParams.c_exc = 0;

dnf_sacParams=initDNFParams();    %Saccades field
dnf_sacParams.tau_u=5;
dnf_sacParams.c_inh = 30;
dnf_sacParams.c_exc = 18;

dnf_rchParams=initDNFParams();    %Reach field
dnf_rchParams.tau_u=5;
dnf_rchParams.c_inh = 30;
dnf_rchParams.c_exc = 18;

NexpRewNeurons = 100;
NcuedNeurons = 100;
Wcued        = zeros(2,NcuedNeurons,dnf_rchParams.fieldSize);
if init_cue_w>0
    % saccade weights
    Wcued(1,1:50,:)=normrnd(1,0.1,50,dnf_sacParams.fieldSize);
    Wcued(1,51:100,:)=max(0,normrnd(0.3,0.1,50,dnf_sacParams.fieldSize));
    % reach weights
    Wcued(2,1:50,:)=max(0,normrnd(.3,0.1,50,dnf_rchParams.fieldSize));
    Wcued(2,51:100,:)=normrnd(1,0.1,50,dnf_rchParams.fieldSize);
else
    Wcued(1,:,:)=max(0,normrnd(0.3,0.1,NcuedNeurons,dnf_sacParams.fieldSize));
    Wcued(2,:,:)=max(0,normrnd(0.3,0.1,NcuedNeurons,dnf_rchParams.fieldSize));
end
for i=1:NcuedNeurons
    for j=1:dnf_rchParams.fieldSize
        Wcued(:,i,j)=Wcued(:,i,j)/sum(Wcued(:,i,j));
    end
end

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Stochastic optimal control framework
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
lqgParams_saccade=initLQGParamsSaccade();
lqgParams_reach  =initLQGParamsReach();

%%%%%%%%%%%%%%%%%%
% Protocols
%%%%%%%%%%%%%%%%%%
protocolParams=initTargetFirstProtocolParams();
protocolParams.target_positions= [12 30];
protocolParams.init_expected_reward = [1.0];

%%%%%%%%%%%%%%%%%%%%%%%%%%
% simulation time course %
%%%%%%%%%%%%%%%%%%%%%%%%%%
simParams=initSimulationParams(protocolParams);
reward_type=0; % 0=none, 1=correct effector, 2=left, 3=right
cue_alpha=0.01;
spatial_alpha=0.005;
lambda1=.01;
lambda2=0.996;
reward_da=1.0;

% set times at which field activities are stored (different variants)
% tStoreFields = [100, 200]; % select specific time steps
tStoreFields = 1:simParams.tMax; % store field activities at every time step
freq=1;

%%%%%%%%%%%%%%%%%%
% initialization %
%%%%%%%%%%%%%%%%%%
dnf_in_e = initDNF(dnf_inputParams, tStoreFields);
dnf_in_h = initDNF(dnf_inputParams, tStoreFields);
dnf_er_e = initDNF(dnf_expRewParams, tStoreFields);
dnf_er_h = initDNF(dnf_expRewParams, tStoreFields);
dnf_sac    = initDNF(dnf_sacParams   , tStoreFields);
dnf_rch    = initDNF(dnf_rchParams   , tStoreFields);


stimulus_input_eye     = zeros(1,dnf_in_e.params.fieldSize);
stimulus_input_hand     = zeros(1,dnf_in_h.params.fieldSize);
stimulus_sac        = zeros(1,dnf_sac.params.fieldSize);
stimulus_rch        = zeros(1,dnf_rch.params.fieldSize);
cue_context         = max(0,normrnd(0.15,0.05,1,NcuedNeurons));
effort_sac          = max(0,normrnd(0.15,0.05,1,dnf_sac.params.fieldSize));
effort_rch          = max(0,normrnd(0.15,0.05,1,dnf_rch.params.fieldSize));
expected_reward = max(0,normrnd(0.15,0.05,NexpRewNeurons,NexpRewNeurons));
expected_reward_body = max(0,normrnd(0.15,0.05,NexpRewNeurons,NexpRewNeurons));
expected_reward_eye = max(0,normrnd(0.15,0.05,1,dnf_sac.params.fieldSize));
expected_reward_hand = max(0,normrnd(0.15,0.05,1,dnf_rch.params.fieldSize));
CostSac =0;
CostRch =0;
elig_reach=zeros(size(dnf_rch.output_u));
elig_saccade=zeros(size(dnf_sac.output_u));
elig_spatial=zeros(NexpRewNeurons,NexpRewNeurons);


%%%%%%%%%%%%%%%%%%%%
% Init Controllers %
%%%%%%%%%%%%%%%%%%%%
protocolTxt=uicontrol(fig,'Style','Text','Units','Norm','Position',[gapWidth, 5*controlFieldHeight, controlFieldWidth/2, controlFieldHeight], 'String', 'Protocol', 'HorizontalAlignment', 'left', 'BackgroundColor', 'w');
protocolDropdown=uicontrol(fig,'Style','Popupmenu', 'Units', 'Norm', 'Position', [gapWidth+.05, 5*controlFieldHeight, controlFieldWidth/2, controlFieldHeight], 'String', 'targets first|cue first|memory', 'Callback', @setProtocol);
cueTxt=uicontrol(fig,'Style','Text','Units','Norm','Position',[controlFieldWidth/2+2*gapWidth+.05, 5*controlFieldHeight, controlFieldWidth/2, controlFieldHeight], 'String', 'Cue', 'HorizontalAlignment', 'left', 'BackgroundColor', 'w');
cueDropdown=uicontrol(fig, 'Style', 'Popupmenu', 'Units', 'Norm', 'Position', [.075+controlFieldWidth/2+3*gapWidth, 5*controlFieldHeight, controlFieldWidth/2, controlFieldHeight], 'String', 'saccade|reach|free choice|saccade and reach', 'Callback', @setCue);

startButton=uicontrol(fig, 'Style', 'Pushbutton', 'String', 'Start', 'Units', 'Norm', 'Position', [.075+controlFieldWidth+4*gapWidth, 5*controlFieldHeight+.0125, controlFieldWidth/4, .025], 'Callback', @startSim);
pauseButton=uicontrol(fig, 'Style', 'Pushbutton', 'String', 'Pause', 'Units', 'Norm', 'Position', [.075+controlFieldWidth+5*gapWidth+controlFieldWidth/4, 5*controlFieldHeight+.0125, controlFieldWidth/4, .025], 'Callback', @pauseSim);
stopButton=uicontrol(fig, 'Style', 'Pushbutton', 'String', 'Stop', 'Units', 'Norm', 'Position', [.075+controlFieldWidth+6*gapWidth+2*controlFieldWidth/4, 5*controlFieldHeight+.0125, controlFieldWidth/4, .025], 'Callback', @stopSim);
saveCueButton=uicontrol(fig, 'Style', 'Pushbutton', 'String', 'Save Cue W', 'Units', 'Norm', 'Position', [.075+controlFieldWidth+7*gapWidth+3*controlFieldWidth/4, 5*controlFieldHeight+.0125, controlFieldWidth/3, .025], 'Callback', @saveCueW);
loadCueButton=uicontrol(fig, 'Style', 'Pushbutton', 'String', 'Load Cue W', 'Units', 'Norm', 'Position', [.075+controlFieldWidth+8*gapWidth+3*controlFieldWidth/4+controlFieldWidth/3, 5*controlFieldHeight+.0125, controlFieldWidth/3, .025], 'Callback', @loadCueW);
saveSpatialButton=uicontrol(fig, 'Style', 'Pushbutton', 'String', 'Save Spatial W', 'Units', 'Norm', 'Position', [.075+controlFieldWidth+7*gapWidth+3*controlFieldWidth/4, 4*controlFieldHeight-gapWidth, controlFieldWidth/3, .025], 'Callback', @saveSpatialW);
loadSpatialButton=uicontrol(fig, 'Style', 'Pushbutton', 'String', 'Load Spatal W', 'Units', 'Norm', 'Position', [.075+controlFieldWidth+8*gapWidth+3*controlFieldWidth/4+controlFieldWidth/3, 4*controlFieldHeight-gapWidth, controlFieldWidth/3, .025], 'Callback', @loadSpatialW);

recordCheckbox=uicontrol(fig, 'Style', 'Checkbox', 'String', 'Record', 'Value', 0, 'Units', 'Norm', 'Position', [.075+controlFieldWidth+4*gapWidth, 4*controlFieldHeight-gapWidth, controlFieldWidth/4, .025], 'Callback', @toggleRecord);
recordFilefield=uicontrol(fig, 'Style', 'Edit', 'String', 'DNF_LQG_simulation.avi', 'Units', 'Norm', 'Position', [.075+controlFieldWidth+5*gapWidth+controlFieldWidth/4, 4*controlFieldHeight-gapWidth, controlFieldWidth/2, .025]);

targetsTxt=uicontrol(fig,'Style','Text','Units','Norm','Position',[gapWidth, 4*controlFieldHeight-gapWidth, controlFieldWidth/2, controlFieldHeight], 'String', 'Targets', 'HorizontalAlignment','left','BackgroundColor','w');
targetsDropdown=uicontrol(fig, 'Style', 'Popupmenu', 'String', '1|2|3', 'Units', 'Norm', 'Position', [gapWidth+.05, 4*controlFieldHeight-gapWidth+.0125, controlFieldWidth/6, .025], 'Callback', @setNumTargets);

trialsTxt=uicontrol(fig,'Style','Text','Units','Norm','Position',[controlFieldWidth/2+2*gapWidth+.05, 4*controlFieldHeight-gapWidth, controlFieldWidth/2, controlFieldHeight], 'String', 'Trials', 'HorizontalAlignment','left','BackgroundColor','w');
trialsDropdown=uicontrol(fig, 'Style', 'Popupmenu', 'String', '1|10|25|50|100|250|500', 'Units', 'Norm', 'Position', [.075+controlFieldWidth/2+3*gapWidth, 4*controlFieldHeight-gapWidth+.0125, controlFieldWidth/6, .025], 'Callback', @setNumTrials);

rewardTxt=uicontrol(fig,'Style','Text','Units','Norm','Position',[gapWidth, 3*controlFieldHeight-2*gapWidth, controlFieldWidth/2, controlFieldHeight], 'String', 'Reward', 'HorizontalAlignment', 'left', 'BackgroundColor', 'w');
rewardDropdown=uicontrol(fig, 'Style', 'Popupmenu', 'String', 'none|correct effector|left target|right target', 'Units', 'Norm', 'Position', [gapWidth+.05, 3*controlFieldHeight-2*gapWidth, controlFieldWidth/2, controlFieldHeight], 'Callback', @setReward);

lesionTxt=uicontrol(fig,'Style','Text','Units','Norm','Position',[gapWidth, 2*controlFieldHeight-3*gapWidth, controlFieldWidth/2, controlFieldHeight], 'String', 'Lesion', 'HorizontalAlignment', 'left', 'BackgroundColor', 'w');
lesionDropdown=uicontrol(fig, 'Style', 'Popupmenu', 'String', 'none|left saccade|right saccade|left reach|right reach', 'Units', 'Norm', 'Position', [gapWidth+.05, 2*controlFieldHeight-3*gapWidth, controlFieldWidth/2, controlFieldHeight], 'Callback', @setLesion);

freqTxt=uicontrol(fig,'Style','Text','Units','Norm','Position',[controlFieldWidth/2+2*gapWidth+.05, 3*controlFieldHeight-2*gapWidth, controlFieldWidth/2, controlFieldHeight], 'String', 'Freq', 'HorizontalAlignment', 'left', 'BackgroundColor','w');
freqDropdown=uicontrol(fig, 'Style', 'Popupmenu', 'String', '1|2|5|10|20|100', 'Units', 'Norm', 'Position', [.075+controlFieldWidth/2+3*gapWidth, 3*controlFieldHeight-2*gapWidth, controlFieldWidth/6, controlFieldHeight], 'Callback', @setFreq);

trialLblTxt=uicontrol(fig,'Style','Text','Units','Norm','Position',[.075+controlFieldWidth+4*gapWidth, 3*controlFieldHeight-2*gapWidth, controlFieldWidth/2, controlFieldHeight], 'String', 'Trial:', 'HorizontalAlignment', 'left', 'BackgroundColor','w');
trialTxt=uicontrol(fig,'Style','Text','Units','Norm','Position',[.1+controlFieldWidth+5*gapWidth, 3*controlFieldHeight-2*gapWidth, controlFieldWidth/2, controlFieldHeight], 'String', '1', 'HorizontalAlignment', 'left', 'BackgroundColor','w');


nControlParams = length(controlParamNames);
sliders = zeros(nControlParams, 1);
textFields = zeros(nControlParams, 1);
if debug>0
    for i = 1 : nControlParams
        eval(['tmp = ' controlParamNames{i} ';']);
        sliders(i) = uicontrol(param_fig, 'Style', 'Slider', 'Units', 'Norm', 'Position', ...
            [controlPosX(i)+textWidth+gapWidth, controlPosY(i), sliderWidth, sliderHeight], ...
            'Value', tmp, 'Min', controlMin(i), 'Max', controlMax(i), 'Callback', @sliderCallback);
        textFields(i) = uicontrol(param_fig,'Style','Text','Units','Norm','HorizontalAlignment', 'left', ...
            'String',[controlParamNames{i} '=' num2str(tmp, textFormat{i})], 'BackgroundColor', 'w',...
            'Position',[controlPosX(i)+gapWidth, controlPosY(i) textWidth controlFieldHeight]);
    end
end


%%%%%%%%%%%%%%
% simulation %
%%%%%%%%%%%%%%
flag_evaluate_sac    = 1;
flag_evaluate_rch    = 1;
cnt_evaluate_sac     = 0;
cnt_evaluate_rch     = 0;
recordSim = 0;

% plot graphs
actPlot_u0_eye=0;
outPlot_u0_eye=0;
actPlot_u0_hand=0;
outPlot_u0_hand=0;
actPlot_u1=0;
outPlot_u1=0;
actPlot_u2=0;
outPlot_u2=0;
inPlot0_eye=0;
inPlot0_hand=0;
inPlot1=0;
inPlot2=0;
kernelPlot0_eye=0;
kernelPlot0_hand=0;
kernelPlot1=0;
kernelPlot2=0;
if debug>0
    
    changeAxes(SpatialInAxes_Sac,fig);
    cla;
    hold on;
    plot([0,dnf_in_e.params.fieldSize-1],[0,0],'Linestyle',':','Linewidth',1);
    actPlot_u0_eye = plot(0:dnf_in_e.params.fieldSize-1,dnf_in_e.field_u,'color','b','Linewidth',3);
    outPlot_u0_eye = plot(0:dnf_in_e.params.fieldSize-1,10*dnf_in_e.output_u,'color','r','Linewidth',1);
    inPlot0_eye = plot(0:dnf_in_e.params.fieldSize-1,stimulus_input_eye+dnf_in_e.params.h_u,'color','g','Linewidth',1);
    set(gca,'ylim',[-15,15],'xlim',[0,dnf_in_e.params.fieldSize-1],'Ytick',[-10,0,10]);
    ylabel('Input Field Eye','Fontsize',12);
    hold off;
    
    changeAxes(SpatialInAxes_Rch,fig);
    cla;
    hold on;
    plot([0,dnf_in_h.params.fieldSize-1],[0,0],'Linestyle',':','Linewidth',1);
    actPlot_u0_hand = plot(0:dnf_in_h.params.fieldSize-1,dnf_in_h.field_u,'color','b','Linewidth',3);
    outPlot_u0_hand = plot(0:dnf_in_h.params.fieldSize-1,10*dnf_in_h.output_u,'color','r','Linewidth',1);
    inPlot0_hand = plot(0:dnf_in_h.params.fieldSize-1,stimulus_input_hand+dnf_in_h.params.h_u,'color','g','Linewidth',1);
    set(gca,'ylim',[-15,15],'xlim',[0,dnf_in_h.params.fieldSize-1],'Ytick',[-10,0,10]);
    ylabel('Input Field Hand','Fontsize',12);    
    
    changeAxes(SacAxes,fig);
    cla;
    hold on;
    plot([0,dnf_sac.params.fieldSize-1],[0,0],'Linestyle',':','Linewidth',1);
    actPlot_u1 = plot(0:dnf_sac.params.fieldSize-1,dnf_sac.field_u,'color','b','Linewidth',3);
    outPlot_u1 = plot(0:dnf_sac.params.fieldSize-1,10*dnf_sac.output_u,'color','r','Linewidth',1);
    inPlot1 = plot(0:dnf_sac.params.fieldSize-1,stimulus_sac+dnf_sac.params.h_u,'color','g','Linewidth',1);
    eligPlot1=plot(0:dnf_sac.params.fieldSize-1,10*elig_saccade,'color','k','Linewidth',1,'LineStyle','--');
    set(gca,'ylim',[-15,15],'xlim',[0,dnf_sac.params.fieldSize-1],'Ytick',[-10,0,10]);
    ylabel('Saccade Field','Fontsize',12);
    hold off;
    
    changeAxes(ReachAxes,fig);
    cla;
    hold on;
    plot([0,dnf_rch.params.fieldSize-1],[0,0],'Linestyle',':','Linewidth',1);
    actPlot_u2 = plot(0:dnf_rch.params.fieldSize-1,dnf_rch.field_u,'color','b','Linewidth',3);
    outPlot_u2 = plot(0:dnf_rch.params.fieldSize-1,10*dnf_rch.output_u,'color','r','Linewidth',1);
    inPlot2 = plot(0:dnf_rch.params.fieldSize-1,stimulus_rch+dnf_rch.params.h_u,'color','g','Linewidth',1);
    eligPlot2=plot(0:dnf_rch.params.fieldSize-1,10*elig_reach,'color','k','Linewidth',1,'LineStyle','--');
    set(gca,'ylim',[-15,15],'xlim',[0,dnf_rch.params.fieldSize-1],'Ytick',[-10,0,10]);
    ylabel('Reach Field','Fontsize',12);
    hold off;
    
    changeAxes(CueAxes,fig);
    cla;
    hold on;
    plot([0,NcuedNeurons-1],[0,0],'Linestyle',':','Linewidth',1);
    outPlot_u3 = plot(0:NcuedNeurons-1,cue_context,'color','r','Linewidth',1);
    set(gca,'ylim',[-.1,1.1],'xlim',[0,NcuedNeurons-1],'Ytick',[-1,0,1]);
    ylabel('Cue Field','Fontsize',12);
    hold off;
    
    
    changeAxes(kernelAxes,fig);
    cla;
    hold on;
    plot([-dnf_sac.params.halfField,dnf_sac.params.halfField],[0,0],'Linestyle',':','Linewidth',1);
    kernelPlot0_eye = plot(-dnf_in_e.params.halfField:dnf_in_e.params.halfField, ...
        [zeros(1, dnf_in_e.params.halfField-dnf_in_e.kSize_uu) dnf_in_e.kernel_uu zeros(1, dnf_in_e.params.halfField-dnf_in_e.kSize_uu)] - dnf_in_e.params.g_inh, 'Color', 'r', 'Linewidth', 3);
    kernelPlot0_hand = plot(-dnf_in_h.params.halfField:dnf_in_h.params.halfField, ...
        [zeros(1, dnf_in_h.params.halfField-dnf_in_h.kSize_uu) dnf_in_h.kernel_uu zeros(1, dnf_in_h.params.halfField-dnf_in_h.kSize_uu)] - dnf_in_h.params.g_inh, 'Color', 'r', 'Linewidth', 3, 'Linestyle', '--');
    kernelPlot1 = plot(-dnf_sac.params.halfField:dnf_sac.params.halfField, ...
        [zeros(1, dnf_sac.params.halfField-dnf_sac.kSize_uu) dnf_sac.kernel_uu zeros(1, dnf_sac.params.halfField-dnf_sac.kSize_uu)] - dnf_sac.params.g_inh, 'Color', 'g', 'Linewidth', 3);
    kernelPlot2 = plot(-dnf_rch.params.halfField:dnf_rch.params.halfField, ...
        [zeros(1, dnf_rch.params.halfField-dnf_rch.kSize_uu) dnf_rch.kernel_uu zeros(1, dnf_rch.params.halfField-dnf_rch.kSize_uu)] - dnf_rch.params.g_inh, 'Color', 'b', 'Linewidth', 3, 'LineStyle', '--');
    set(gca,'ylim',[-2,2],'xlim',[-dnf_rch.params.halfField,dnf_rch.params.halfField],'Ytick',[-2,-1,0,1,2]);
    ylabel('Interaction Kernel','Fontsize',12);
    hold off;
    
    changeAxes(EffortSacAxes,fig);
    cla;
    hold on;
    plot([0,dnf_sac.params.fieldSize-1],[0,0],'Linestyle',':','Linewidth',1);
    outPlot_u4 = plot(0:dnf_sac.params.fieldSize-1,effort_sac,'color','r','Linewidth',1);
    set(gca,'ylim',[-.2,1.5],'xlim',[0,dnf_sac.params.fieldSize-1],'Ytick',[-1,0,1]);
    ylabel('Effort Saccade','Fontsize',12);
    hold off;
    
    changeAxes(EffortRchAxes,fig);
    cla;
    hold on;
    plot([0,dnf_rch.params.fieldSize-1],[0,0],'Linestyle',':','Linewidth',1);
    outPlot_u5 = plot(0:dnf_rch.params.fieldSize-1,effort_sac,'color','r','Linewidth',1);
    set(gca,'ylim',[-.2,1.5],'xlim',[0,dnf_rch.params.fieldSize-1],'Ytick',[-1,0,1]);
    ylabel('Effort Reach','Fontsize',12);
    hold off;
    
    changeAxes(CueWAxes,fig);
    cla;
    hold on;
    plot([0,dnf_sac.params.fieldSize-1],[0,0],'Linestyle',':','Linewidth',1);
    outPlot_w1 = plot(0:dnf_sac.params.fieldSize-1,squeeze(mean(Wcued(1,1:50,:))),'color','r','Linewidth',1);
    outPlot_w2 = plot(0:dnf_sac.params.fieldSize-1,squeeze(mean(Wcued(1,51:100,:))),'color','g','Linewidth',1);
    outPlot_w3 = plot(0:dnf_rch.params.fieldSize-1,squeeze(mean(Wcued(2,1:50,:))),'color','r','Linewidth',1,'LineStyle','--');
    outPlot_w4 = plot(0:dnf_rch.params.fieldSize-1,squeeze(mean(Wcued(2,51:100,:))),'color','g','Linewidth',1,'LineStyle','--');
    set(gca, 'ylim',[0,1],'xlim',[0,dnf_rch.params.fieldSize-1]);
    ylabel('Cue W','Fontsize',12);
    hold off;
    
    changeAxes(RewardEgoAxes,fig);
    cla;
    hold on;
    plot([0,dnf_rch.params.fieldSize-1],[0,0],'Linestyle',':','Linewidth',1);
    outPlot_u8 = plot(0:dnf_er_e.params.fieldSize-1,dnf_er_e.output_u,'color','r','Linewidth',1);
    outPlot_u9 = plot(0:dnf_er_h.params.fieldSize-1,dnf_er_h.output_u,'color','g','Linewidth',1);
    set(gca,'ylim',[0 1.1],'xlim',[0,dnf_er_e.params.fieldSize-1],'Ytick',[0,1]);
    ylabel('E(Reward) - Effector','Fontsize',12);
    hold off;

    changeAxes(RewardAlloAxes,fig);
    outPlot_er1=imagesc(expected_reward');
    set(gca,'ylim',[1,NexpRewNeurons],'xlim',[1,NexpRewNeurons],'Ydir','normal','Xtick',[20 60 100],'clim',[0 1]);
    cb=colorbar();
    set(cb,'Ytick',[0,1]);
    ylabel('E(Reward) - Body');

    changeAxes(RewardBodyAxes,fig);
    outPlot_er2=imagesc(expected_reward_body');
    set(gca,'ylim',[1,NexpRewNeurons],'xlim',[1,NexpRewNeurons],'Ydir','normal','Xtick',[20 60 100],'clim',[0 1]);
    cb=colorbar();
    set(cb,'Ytick',[-.1,1.1]);
end

cnt_movie = 0;

trial=1;

expected_reward=initExpectedReward(NexpRewNeurons,protocolParams);

while true
    if simParams.running>0
        stim_input_eye=zeros(1,dnf_in_e.params.fieldSize);
        stim_input_hand=zeros(1,dnf_in_h.params.fieldSize);
        
        history=initHistory(simParams.nTrials, freq, simParams.tMax, dnf_in_e, dnf_in_h, dnf_sac, dnf_rch,NcuedNeurons,NexpRewNeurons);                

        % loop over trials
        while trial<=simParams.nTrials
            while simParams.running<1
                pause(.001);
            end
            plotted=0;
            disp(['trial ' num2str(trial)]);
            
            if reward_type==1
                initNumTargets(1);
            end
            protocolParams=resetProtocolParams(protocolParams);
            simParams=resetSimParams(simParams, protocolParams);
            flag_evaluate_sac=1;
            flag_evaluate_rch=1;
            cnt_evaluate_sac= 0;
            cnt_evaluate_rch= 0;

            dnf_in_e = resetDNF(dnf_in_e, tStoreFields);
            dnf_in_h = resetDNF(dnf_in_h, tStoreFields);
            dnf_er_e = resetDNF(dnf_in_e, tStoreFields);
            dnf_er_h = resetDNF(dnf_in_h, tStoreFields);
            dnf_sac    = resetDNF(dnf_sac, tStoreFields);
            dnf_rch    = resetDNF(dnf_rch, tStoreFields);
            
            SaccadicMovement=zeros(2,simParams.tMax);
            ReachingMovement=zeros(2,simParams.tMax);
            changeAxes(BehaviorAxes,fig);
            
            dopamine_reward=0;
            elig_reach=zeros(size(dnf_rch.output_u));
            elig_saccade=zeros(size(dnf_sac.output_u));
            elig_spatial=zeros(NexpRewNeurons,NexpRewNeurons);
            
            if protocolParams.random_cue>0
                if rand()>=.5
                    protocolParams.cue=0;
                else
                    protocolParams.cue=1;
                end
            end
            
            %Initialization of the saccade and reach state
            SaccadicStateTotal = [simParams.Xeye(1:4) zeros(4,size(lqgParams_saccade.Q,3))];
            ReachingStateTotal = [simParams.Xhand(1:8) zeros(8,size(lqgParams_reach.Q,3))];
            ControlSacX = 0; ControlSacY = 0; ControlRchX = 0; ControlRchY = 0;
            
            trial_cue=zeros(simParams.tMax,NcuedNeurons);
            trial_effort_sac=zeros(simParams.tMax,dnf_sac.params.fieldSize);
            trial_effort_rch=zeros(simParams.tMax,dnf_rch.params.fieldSize);

            % loop over time steps
            for t = 1 : simParams.tMax
                while simParams.running<1
                    pause(.001);
                end
                if simParams.flag_contact < 1
                    if rem(t,simParams.revaluate_threshold) == 0  % Re-evaluate the reaching/saccade policies
                        cnt_evaluate_sac = 0;
                        cnt_evaluate_rch = 0;
                        [MotorOutSac CostSac ActivePopSac] = LQGdnf(dnf_sac,simParams.output_u_sac_thr,...
                            lqgParams_saccade,simParams.Xeye,simParams.minEyeDist, ...
                            protocolParams.target_positions,0);
                        if isempty(ActivePopSac) == 0
                            [SaccadicStateTotal ControlSacX ControlSacY] = ActionSelection(MotorOutSac,ActivePopSac,dnf_sac,lqgParams_saccade,simParams.Xeye,0);
                            ControlSacX_previous = ControlSacX;
                            ControlSacY_previous = ControlSacY;
                        elseif isempty(ActivePopSac) == 0 & sum(ControlSacX_previous)>0 % If there is no active population while moving, continue the same action
                            SaccadicStateTotal = RunSimTraj(lqgParams_saccade,[ControlSacX_previous(simParams.revaluate_threshold)*ones(size(lqgParams_saccade.Q,3)-1,1) ControlSacY_previous(simParams.revaluate_threshold)*ones(size(lqgParams_saccade.Q,3)-1,1)]',simParams.Xeye,0);
                        end
                            
                        [MotorOutRch CostRch ActivePopRch] = LQGdnf(dnf_rch,simParams.output_u_rch_thr,...
                            lqgParams_reach,simParams.Xhand,simParams.minHandDist, ...
                            protocolParams.target_positions,1);
                        if isempty(ActivePopRch) == 0
                            [ReachingStateTotal ControlRchX ControlRchY]  = ActionSelection(MotorOutRch,ActivePopRch,dnf_rch,lqgParams_reach,simParams.Xhand,1);
                            ControlRchX_previous = ControlRchX;
                            ControlRchY_previous = ControlRchY;
                        elseif isempty(ActivePopRch) == 0 & sum(ControlRchX_previous)>0 % If there is no active population while moving, continue the same action
                            ReachingStateTotal = RunSimTraj(lqgParams_reach,[ControlRchX_previous(simParams.revaluate_threshold)*ones(size(lqgParams_reach.Q,3)-1,1) ControlRchY_previous(simParams.revaluate_threshold)*ones(size(lqgParams_reach.Q,3)-1,1)]',simParams.Xhand,1);
                        end  
                    end
                    
                    if max(dnf_sac.output_u) > simParams.thr_perf_saccade & isempty(ActivePopSac) ==0  % Start the saccadic movements
                        if abs(sum(SaccadicMovement(:)))==0 && cnt_evaluate_sac == 0
                            disp('starting saccade');
                            simParams.chosen_effector=1;
                        end
                        cnt_evaluate_sac = cnt_evaluate_sac + 1;
                        SaccadicMovement(:,t) = SaccadicStateTotal([1,3],cnt_evaluate_sac);
                        simParams.Xeye(1:end-2) = SaccadicStateTotal(:,cnt_evaluate_sac);
                        simParams.xe   = SaccadicMovement(1,t);
                        simParams.ye   = SaccadicMovement(2,t);
                        simParams.eye_dist_targets=[];
                        for i=1:size(protocolParams.target_positions,1)
                            rx=protocolParams.target_positions(i,1);
                            ry=protocolParams.target_positions(i,2);
                            simParams.eye_dist_targets=[simParams.eye_dist_targets; sqrt((SaccadicMovement(1,t)-rx)^2  + (SaccadicMovement(2,t)-ry)^2)];
                        end                
                        simParams.minEyeDist=min(simParams.eye_dist_targets);
                        changeAxes(BehaviorAxes,fig);
                        hold on;
                        plot(SaccadicMovement(1,t),SaccadicMovement(2,t),'o','color',[1 0 0]);
                    end
                    
                    
                    if max(dnf_rch.output_u) > simParams.thr_perf_reach & isempty(ActivePopRch) ==0
                        if abs(sum(ReachingMovement(:)))==0 && cnt_evaluate_rch == 0
                            disp('starting reach');
                            simParams.chosen_effector=2;
                        end
                        cnt_evaluate_rch = cnt_evaluate_rch + 1;
                        ReachingMovement(:,t) = ReachingStateTotal([1,2],cnt_evaluate_rch);
                        simParams.Xhand(1:end-2) = ReachingStateTotal(:,cnt_evaluate_rch);
                        simParams.xh    = ReachingMovement(1,t);
                        simParams.yh    = ReachingMovement(2,t);
                        simParams.hand_dist_targets=[];
                        for i=1:size(protocolParams.target_positions,1)
                            rx=protocolParams.target_positions(i,1);
                            ry=protocolParams.target_positions(i,2);
                            simParams.hand_dist_targets=[simParams.hand_dist_targets; sqrt((ReachingMovement(1,t)-rx)^2  + (ReachingMovement(2,t)-ry)^2)];
                        end
                        simParams.minHandDist=min(simParams.hand_dist_targets);
                        changeAxes(BehaviorAxes,fig);
                        hold on;
                        plot(ReachingMovement(1,t),ReachingMovement(2,t),'o','color',[0 1 0]);
                    end
                else
                    CostSac = 0;
                    CostRch = 0;
                end

                chosen_target_idx=0;
                if simParams.minEyeDist < simParams.Threshold_contact
                    simParams.flag_contact = 1;
                    chosen_target_idx=find(simParams.eye_dist_targets==simParams.minEyeDist);
                    if simParams.xe<-1
                        simParams.chosen_target=1;
                    elseif simParams.xe>1
                        simParams.chosen_target=2;
                    end
                elseif simParams.minHandDist < simParams.Threshold_contact
                    simParams.flag_contact = 1;
                    chosen_target_idx=find(simParams.hand_dist_targets==simParams.minHandDist);
                    if simParams.xh<-1
                        simParams.chosen_target=1;
                    elseif simParams.xh>1
                        simParams.chosen_target=2;
                    end
                end
                
		changeAxes(BehaviorAxes,fig)
		hold on;

                if protocolParams.protocol == 0 % target first
                    [protocolParams stim_input_eye stim_input_hand cue_context]=runTargetFirstProtocol(t, simParams,...
                       protocolParams, dnf_in_e, dnf_in_h);
                elseif protocolParams.protocol == 1 % cue first
                    [protocolParams stim_input_eye stim_input_hand cue_context]=runCueFirstProtocol(t, simParams,...
                    protocolParams, dnf_in_e, dnf_in_h);
                elseif protocolParams.protocol == 2 % memory
                    [protocolParams stim_input_eye stim_input_hand cue_context]=runMemoryProtocol(t, simParams,...
                    protocolParams, dnf_in_e, dnf_in_h);
                end

                hold off;
                
                disp(['t=' num2str(t)]);
                
                stimulus_field_eye     = max(0,normrnd(0.15,0.05,1,dnf_in_e.params.fieldSize))+2*stim_input_eye;
                stimulus_field_hand    = max(0,normrnd(0.15,0.05,1,dnf_in_h.params.fieldSize))+2*stim_input_hand;
                effort_sac      = normrnd(0,0.05,1,dnf_sac.params.fieldSize) + 10*CostSac;
                effort_rch      = normrnd(0,0.05,1,dnf_rch.params.fieldSize) + 10*CostRch;
                
                expected_reward_body=normrnd(0,0.05,NexpRewNeurons,NexpRewNeurons);
                expected_reward_hand=normrnd(0,0.05,1,dnf_in_h.params.fieldSize);
                expected_reward_eye=normrnd(0,0.05,1,dnf_in_h.params.fieldSize);
                if max(stim_input_eye) > 0 || max(stim_input_hand) > 0
                    expected_reward_body=expected_reward_body+expected_reward; 
                end
                expected_reward_eye=expected_reward_eye+convertTwoDReferenceFramePolar(expected_reward_body, ...
                    dnf_in_e.params.fieldSize, NexpRewNeurons, [simParams.xe simParams.ye],...
                    protocolParams.target_positions);
                expected_reward_hand=expected_reward_hand+convertTwoDReferenceFramePolar(expected_reward_body, ...
                    dnf_in_h.params.fieldSize, NexpRewNeurons, [simParams.xh simParams.yh],...
                    protocolParams.target_positions);
                
                inhib_w=-1.25;
                context_w=0.1;
                in_w=8.0;
                effort_w=-0.1;
                er_w=0.2;
                lesion_w=-0.1;

                lesion_sac=zeros(1,dnf_in_e.params.fieldSize);
                lesion_rch=zeros(1,dnf_in_h.params.fieldSize);
                if protocolParams.lesion==1 % left eye
                    lesion_sac(1:round(dnf_in_e.params.fieldSize/2))=lesion_w;
                elseif protocolParams.lesion==2 % right eye
                    lesion_sac(round(dnf_in_e.params.fieldSize/2):dnf_in_e.params.fieldSize)=lesion_w;
                elseif protocolParams.lesion==3 % left hand
                    lesion_rch(1:round(dnf_in_h.params.fieldSize/2))=lesion_w;
                elseif protocolParams.lesion==4 % right hand
                    lesion_rch(round(dnf_in_h.params.fieldSize/2):dnf_in_h.params.fieldSize)=lesion_w;
                end

		% input to the saccade motor plan formation DNF
                stimulus_sac = inhib_w*sum(dnf_rch.output_u(:))*ones(1,dnf_rch.params.fieldSize)+...
                    context_w*cue_context*squeeze(Wcued(1,:,:))+in_w*dnf_in_e.output_u+...
                    effort_w*effort_sac+er_w*dnf_er_e.output_u+lesion_sac;

		% input to the reach motor plan formation DNF
                stimulus_rch = inhib_w*sum(dnf_sac.output_u(:))*ones(1,dnf_sac.params.fieldSize)+...
                   context_w*cue_context*squeeze(Wcued(2,:,:))+in_w*dnf_in_h.output_u+...
                   effort_w*effort_rch+er_w*dnf_er_h.output_u+lesion_rch;
                
		% run eye stimulus input field
                dnf_in_e=runDNF(dnf_in_e, stimulus_field_eye, tStoreFields, t);
		% run hand stimulus input field
                dnf_in_h=runDNF(dnf_in_h, stimulus_field_hand, tStoreFields, t);
		% run eye expected reward field
                dnf_er_e=runDNF(dnf_er_e, 10*expected_reward_eye, tStoreFields, t);
		% run hand expected reward field
                dnf_er_h=runDNF(dnf_er_h, 10*expected_reward_hand, tStoreFields, t);
		% run saccade motor plan formation DNF
                dnf_sac=runDNF(dnf_sac, stimulus_sac, tStoreFields, t);
		% run reach motor plan formation DNF
                dnf_rch=runDNF(dnf_rch, stimulus_rch, tStoreFields, t);
                
		% compute reach eligibility (decaying copy of reach motor plan formation DNF output)
                elig_reach=min(1,elig_reach+lambda1*dnf_rch.output_u-(1-lambda2)*elig_reach);
		% compute saccade eligibility (decaying copy of reach motor plan formation DNF output)
                elig_saccade=min(1,elig_saccade+lambda1*dnf_sac.output_u-(1-lambda2)*elig_saccade);
                
                trial_cue(t,:)=cue_context;
                trial_effort_sac(t,:)=effort_sac;
                trial_effort_rch(t,:)=effort_rch;

                if debug>0
                    set(inPlot0_eye, 'Ydata', stimulus_field_eye+dnf_in_e.params.h_u);
                    set(inPlot0_hand, 'Ydata', stimulus_field_hand+dnf_in_h.params.h_u);
                    set(actPlot_u0_eye,'Ydata',dnf_in_e.field_u);
                    set(actPlot_u0_hand,'Ydata',dnf_in_h.field_u);
                    set(outPlot_u0_eye,'Ydata',10*dnf_in_e.output_u);
                    set(outPlot_u0_hand,'Ydata',10*dnf_in_h.output_u);
                    set(inPlot1, 'Ydata', stimulus_sac+dnf_sac.params.h_u);
                    set(eligPlot1, 'Ydata', 10*elig_saccade);
                    set(actPlot_u1,'Ydata',dnf_sac.field_u);
                    set(outPlot_u1,'Ydata',10*dnf_sac.output_u);
                    set(inPlot2, 'Ydata', stimulus_rch+dnf_rch.params.h_u);
                    set(eligPlot2, 'Ydata', 10*elig_reach);
                    set(actPlot_u2,'Ydata',dnf_rch.field_u);
                    set(outPlot_u2,'Ydata',10*dnf_rch.output_u);
                    set(outPlot_u3,'Ydata',cue_context);
                    set(outPlot_u4,'Ydata',effort_sac);
                    set(outPlot_u5,'Ydata',effort_rch);
                    set(outPlot_u8,'Ydata',dnf_er_e.output_u);
                    set(outPlot_u9,'Ydata',dnf_er_h.output_u);
                    set(outPlot_er1,'Cdata',expected_reward');
                    set(outPlot_er2,'Cdata',expected_reward_body');
                    drawnow;
                end
                
                if recordSim>0
                    ff(t) = getframe(fig);
                end
            end
            
            % spatial eligibility field (multivariate Gaussian of location moved to with either effector)
            elig_spatial=zeros(NexpRewNeurons,NexpRewNeurons);

            % if an effector contacted a target on this trial
            if simParams.flag_contact>0

		% compute spatial eligibility field
                for x=1:NexpRewNeurons
                    for y=1:NexpRewNeurons
                        % x coordinate for index x, eligspatial is 100x100 and codes for x=-45 to x=45
                        x_coord=(x-NexpRewNeurons/2.0+1.0)/(NexpRewNeurons)*(45.0--45.0);
                        % y coordinate for index y, eligspatial is 100x100 and codes for y=0 to x=50
                        y_coord=(y-1.0)*(50-0.0)/(NexpRewNeurons-1);
			% the eligibility trace is a mutivariate Guassian of the x and y coordinate, centered on the coordinate of the chosen target
                        elig_spatial(x,y)=gauss(x_coord,protocolParams.target_positions(chosen_target_idx,1),5)*gauss(y_coord,protocolParams.target_positions(chosen_target_idx,2),5);
                    end
                end

		% determine if reward was earned
                if (reward_type==1 && ((simParams.chosen_effector==1 && protocolParams.cue==0) || (simParams.chosen_effector==2 && protocolParams.cue==1))) || (reward_type==2 && simParams.chosen_target==1) || (reward_type==3 && simParams.chosen_target==2)
                    dopamine_reward=reward_da;
                end
            end

            % Update expected reward and normalize
            expected_reward=expected_reward+elig_spatial.*(spatial_alpha*dopamine_reward);
            expected_reward=expected_reward./max(expected_reward(:));
            
            % Update and normalize cue weights
            for i=1:NcuedNeurons
                Wcued(1,i,:)=squeeze(Wcued(1,i,:))'+cue_alpha*dopamine_reward*cue_context(i)*elig_saccade;
                Wcued(2,i,:)=squeeze(Wcued(2,i,:))'+cue_alpha*dopamine_reward*cue_context(i)*elig_reach;
                for j=1:dnf_rch.params.fieldSize
                    Wcued(:,i,j)=Wcued(:,i,j)/sum(Wcued(:,i,j));
                end
            end
            
            history=recordHistory(history, trial, ReachingMovement', SaccadicMovement', dnf_in_e, dnf_in_h, dnf_sac, dnf_rch, trial_cue, trial_effort_sac, trial_effort_rch, dopamine_reward, Wcued, expected_reward);

            changeAxes(BehaviorAxes,fig);
            hold off;
            plot(0,0,'.');
            hold on;
            col_incr=1.0/(trial+1);
            dot_width=.5;
            dot_height=.5;
            for j=1:length(history.record_trials)
                trial_idx=history.record_trials(j);
                if trial_idx<=trial
                    symbol='.';
                    if history.dopamine(trial_idx)>0
                        symbol='o';
                    end
                    back_level=1-trial_idx*col_incr;
                    plot(squeeze(history.reaching_movements(j,:,1)), squeeze(history.reaching_movements(j,:,2)), symbol, 'color', [back_level 1 back_level]);
                    plot(squeeze(history.saccade_movements(j,:,1)), squeeze(history.saccade_movements(j,:,2)), symbol, 'color', [1 back_level back_level]);
                end
            end
            hold off;
            set(gca,'ylim',[0,50],'xlim',[-45,45]);
            daspect([1 1 1])
            
            if debug>0
                set(outPlot_w1,'Ydata',squeeze(mean(Wcued(1,1:50,:))));
                set(outPlot_w2,'Ydata',squeeze(mean(Wcued(1,51:100,:))));
                set(outPlot_w3,'Ydata',squeeze(mean(Wcued(2,1:50,:))));
                set(outPlot_w4,'Ydata',squeeze(mean(Wcued(2,51:100,:))));
                drawnow;
            end
            trial=trial+1;
            set(trialTxt, 'String', num2str(trial));
        end
        if plotted<1
            cue='r';
            if protocolParams.cue>0
                cue='g';
            end
            plotHistory(history,500,protocolParams.target_positions,cue,[]);
            plotted=1;
        end
        if recordSim>0 && exist('ff')
            disp('Writing movie file....');
            recordFile=get(recordFilefield,'String');
            movie2avi(ff, recordFile, 'compression', 'None');
            disp('done!');
            clear ff;
        end
        break
    end
    pause(0.001);
end
save('history','history');
disp('Done')
% nStoredFields = simParams.nTrials * length(tStoreFields);
%
% % view evolution of field activities in each trial as mesh plot
% nFieldsPerTrial = length(tStoreFields);
% if 0
%     disp('Press any key to iterate through trials');
%     figure;
%     for i = 1 : nStoredFields
%         plot(0:dnf_sac.params.fieldSize-1, zeros(1, dnf_sac.params.fieldSize), ':k', ...
%             1:dnf_sac.params.fieldSize, dnf_sac.history_s(i, :), '--g', ...
%             1:dnf_sac.params.fieldSize, dnf_sac.history_mu(i, :) + dnf_sac.params.h_u, '-c', ...
%             1:dnf_sac.params.fieldSize, dnf_sac.history_u(i, :), '-b');
%         set(gca, 'XLim', [0 dnf_sac.params.fieldSize-1], 'YLim', [-15 15]);
%         ylabel('activity u');
%         drawnow;
%         pause(0.01);
%     end
%     for i = 1 : nTrials
%         subplot(2, 1, 1);
%         mesh(1:dnf_sac.params.fieldSize, tStoreFields, ...
%             dnf_sac.history_u((i-1)*nFieldsPerTrial+1 : i*nFieldsPerTrial, :));
%         zlabel('activity u');
%         subplot(2, 1, 2);
%         mesh(1:dnf_sac.params.fieldSize, tStoreFields, ...
%             dnf_sac.history_mu((i-1)*nFieldsPerTrial+1 : i*nFieldsPerTrial, :));
%         zlabel('memory trace u');
%         pause
%     end
% end
%
% % view mesh plot of all stored field activities together
% if 0
%     figure;
%     subplot(2, 1, 1);
%     mesh(1:dnf_sac.params.fieldSize, 1:nStoredFields, dnf_sac.history_u(:, :));
%     zlabel('activity u');
%     subplot(2, 1, 2);
%     mesh(1:dnf_sac.params.fieldSize, 1:nStoredFields, dnf_sac.history_mu(:, :));
%     zlabel('memory trace u');
% end

%update paramter values after slider changed
    function sliderCallback(hObject, eventdata) %#ok<INUSD>
        sliderChanged = find(hObject == sliders);
        
        paramName = controlParamNames{sliderChanged};
        tmp = get(sliders(sliderChanged), 'Value');
        set(textFields(sliderChanged), 'String', [paramName '=' num2str(tmp, textFormat{sliderChanged})]);
        eval([paramName '= tmp;']);
        dnf_in_e.kernel_uu = dnf_in_e.params.c_exc * gaussNorm(-dnf_in_e.params.halfField:dnf_in_e.params.halfField, 0, dnf_in_e.params.sigma_exc) ...
            - dnf_in_e.params.c_inh * gaussNorm(-dnf_in_e.params.halfField:dnf_in_e.params.halfField, 0, dnf_in_e.params.sigma_inh) - dnf_in_e.params.g_inh;
        dnf_in_e.kernel_mu = dnf_in_e.params.c_mu * gaussNorm(-dnf_in_e.params.halfField:dnf_in_e.params.halfField, 0, dnf_in_e.params.sigma_mu);
        set(kernelPlot0_eye,'YData',[zeros(1, dnf_in_e.params.halfField-dnf_in_e.kSize_uu) dnf_in_e.kernel_uu zeros(1, dnf_in_e.params.halfField-dnf_in_e.kSize_uu)] - dnf_in_e.params.g_inh);
        dnf_in_h.kernel_uu = dnf_in_h.params.c_exc * gaussNorm(-dnf_in_h.params.halfField:dnf_in_h.params.halfField, 0, dnf_in_h.params.sigma_exc) ...
            - dnf_in_h.params.c_inh * gaussNorm(-dnf_in_h.params.halfField:dnf_in_h.params.halfField, 0, dnf_in_h.params.sigma_inh) - dnf_in_h.params.g_inh;
        dnf_in_h.kernel_mu = dnf_in_h.params.c_mu * gaussNorm(-dnf_in_h.params.halfField:dnf_in_h.params.halfField, 0, dnf_in_h.params.sigma_mu);
        set(kernelPlot0_hand,'YData',[zeros(1, dnf_in_h.params.halfField-dnf_in_h.kSize_uu) dnf_in_h.kernel_uu zeros(1, dnf_in_h.params.halfField-dnf_in_h.kSize_uu)] - dnf_in_h.params.g_inh);
        dnf_sac.kernel_uu = dnf_sac.params.c_exc * gaussNorm(-dnf_sac.params.halfField:dnf_sac.params.halfField, 0, dnf_sac.params.sigma_exc) ...
            - dnf_sac.params.c_inh * gaussNorm(-dnf_sac.params.halfField:dnf_sac.params.halfField, 0, dnf_sac.params.sigma_inh) - dnf_sac.params.g_inh;
        dnf_sac.kernel_mu = dnf_sac.params.c_mu * gaussNorm(-dnf_sac.params.halfField:dnf_sac.params.halfField, 0, dnf_sac.params.sigma_mu);
        set(kernelPlot1,'YData',[zeros(1, dnf_sac.params.halfField-dnf_sac.kSize_uu) dnf_sac.kernel_uu zeros(1, dnf_sac.params.halfField-dnf_sac.kSize_uu)] - dnf_sac.params.g_inh);
        dnf_rch.kernel_uu = dnf_rch.params.c_exc * gaussNorm(-dnf_rch.params.halfField:dnf_rch.params.halfField, 0, dnf_rch.params.sigma_exc) ...
            - dnf_rch.params.c_inh * gaussNorm(-dnf_rch.params.halfField:dnf_rch.params.halfField, 0, dnf_rch.params.sigma_inh) - dnf_rch.params.g_inh;
        dnf_rch.kernel_mu = dnf_rch.params.c_mu * gaussNorm(-dnf_rch.params.halfField:dnf_rch.params.halfField, 0, dnf_rch.params.sigma_mu);
        set(kernelPlot2,'YData',[zeros(1, dnf_rch.params.halfField-dnf_rch.kSize_uu) dnf_rch.kernel_uu zeros(1, dnf_rch.params.halfField-dnf_rch.kSize_uu)] - dnf_rch.params.g_inh);
    end

    function setLesion(hObject, eventdata)
        protocolParams.lesion=get(hObject, 'Value')-1;
    end

    function setProtocol(hObject, eventdata)
        old_cue=protocolParams.cue;
        old_positions=protocolParams.target_positions;
        old_exp_reward=protocolParams.init_expected_reward;
        old_lesion=protocolParams.lesion;
        protocol=get(hObject, 'Value')-1;
        if protocol==0
            protocolParams=initTargetFirstProtocolParams();
        elseif protocol==1
            protocolParams=initCueFirstProtocolParams();
        elseif protocol==2
            protocolParams=initMemoryProtocolParams();
        end
        protocolParams.protocol=protocol;
        protocolParams.cue=old_cue;
        protocolParams.target_positions=old_positions;
        protocolParams.init_expected_reward=old_exp_reward;
        protocolParams.lesion=old_lesion;
        oldNTrials=simParams.nTrials;
        simParams=initSimulationParams(protocolParams);
        simParams.nTrials=oldNTrials;
    end

    function initNumTargets(num_targets)
        if num_targets==1
            p=rand();
            if p<.33
                protocolParams.target_positions= [12 30];
            elseif p<.67
                protocolParams.target_positions= [0 30];
            else
                protocolParams.target_positions= [-12 30];
            end
            protocolParams.init_expected_reward = [1.0];
        elseif num_targets==2
            protocolParams.target_positions=[[-12 30]; [12 30]];
            protocolParams.init_expected_reward = [0.5; 0.5];
        elseif num_targets==3
            protocolParams.target_positions=[[-12 30]; [0 30]; [12 30]];
            protocolParams.init_expected_reward = [0.33; 0.33; 0.33];
        end
        oldNTrials=simParams.nTrials;
        oldRunning=simParams.running;
        simParams=initSimulationParams(protocolParams);
        simParams.nTrials=oldNTrials;
        simParams.running=oldRunning;
    end

    function setNumTargets(hObject, eventdata)
        num_targets=get(hObject, 'Value');
        initNumTargets(num_targets);
        expected_reward=initExpectedReward(NexpRewNeurons,protocolParams);
    end

    function setNumTrials(hObject, eventdata)
        num_trials=get(hObject,'Value');
        if num_trials==1
            simParams.nTrials=1;
        elseif num_trials==2
            simParams.nTrials=10;
        elseif num_trials==3
            simParams.nTrials=25;
        elseif num_trials==4
            simParams.nTrials=50;
        elseif num_trials==5
            simParams.nTrials=100;
        elseif num_trials==6
            simParams.nTrials=250;
        elseif num_trials==7
            simParams.nTrials=500;
        end
    end

    function setFreq(hObject, eventdata)
        frequency=get(hObject,'Value');
        if frequency==1
            freq=1;
        elseif frequency==2
            freq=2;
        elseif frequency==3
            freq=5;
        elseif frequency==4
            freq=10;
        elseif frequency==5
            freq=20;
        elseif frequency==6
            freq=100;
        end
    end

    function setReward(hObject, eventdata)
        reward_type=get(hObject, 'Value')-1;
        if reward_type==1
            protocolParams.target_positions=[[-12 30]; [0 30]; [12 30]];
            protocolParams.init_expected_reward = [0.33; 0.33; 0.33];
        end
    end

    function setCue(hObject, eventdata)
        cueName=get(hObject, 'Value');
        if cueName<4
            protocolParams.cue=cueName-1;
            protocolParams.random_cue=0;
        elseif cueName==4
            protocolParams.random_cue=1;
        end
    end

    function saveCueW(hObject, eventdata)
        save('cue_w', 'Wcued');
    end

    function loadCueW(hObject, eventdata)
        A=load('cue_w');
        Wcued=A.Wcued;
    end

    function saveSpatialW(hObject, eventdata)
        save('spatial_w', 'expected_reward');
    end

    function loadSpatialW(hObject, eventdata)
        A=load('spatial_w');
        expected_reward=A.expected_reward;
    end

    function startSim(hObject, eventdata)
        simParams.running=1;
    end

    function pauseSim(hObject, eventdata)
        simParams.running=0;
    end

    function stopSim(hObject, eventdata)
        protocolParams=resetProtocolParams(protocolParams);
        simParams=resetSimParams(simParams, protocolParams);
        flag_evaluate_sac=1;
        flag_evaluate_rch=1;
        cnt_evaluate_sac= 0;
        cnt_evaluate_rch= 0;
        dnf_in_e = resetDNF(dnf_in_e,tStoreFields);
        dnf_in_h = resetDNF(dnf_in_h,tStoreFields);
        dnf_sac    = resetDNF(dnf_sac,tStoreFields);
        dnf_rch    = resetDNF(dnf_rch,tStoreFields);
        clear SaccadicMovement;
        clear ReachingMovement;
        trial=1;
        t=1;
        simParams.running=0;
    end

    function toggleRecord(hObject, eventdata)
        recordValue=get(hObject,'Value');
        recordSim=recordValue;
    end

    function changeFigure(h)
        set(0,'CurrentFigure',h);
    end

    function changeAxes(h,fig)
        set(fig,'CurrentAxes',h);
    end

end



