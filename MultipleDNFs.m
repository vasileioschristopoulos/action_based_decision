function MultipleDNFs(debug)

close all
clc

fig = figure('Position',[50,50,1200,700],'Name','Neural Field',...
    'Color','w','NumberTitle','off','MenuBar','none');


% create axes for field plots
SpatialInAxes = axes('Position',[0.05 0.88 0.40 0.1]);
SacAxes       = axes('Position',[0.05 0.73 0.40 0.1]);
ReachAxes     = axes('Position',[0.05 0.58 0.40 0.1]);
CueAxes       = axes('Position',[0.05 0.43 0.40 0.1]);
kernelAxes    = axes('Position',[0.05 0.28 0.40 0.1]);
EffortSacAxes = axes('Position',[0.50 0.88 0.20 0.1]);
EffortRchAxes = axes('Position',[0.75 0.88 0.20 0.1]);
BehaviorAxes  = axes('Position',[0.50 0.28 0.45 0.4]);


% %create sliders for model parameters
controlFieldHeight = 0.04;
controlFieldWidth = 0.3;
sliderWidth = 0.15;
gapWidth = 0.01;
textWidth = controlFieldWidth - sliderWidth - gapWidth;

controlParamNames = {'dnf_sac.params.h_u','dnf_sac.params.c_exc','dnf_sac.params.c_inh','dnf_sac.params.g_inh','dnf_sac.params.q_u',...
    'dnf_sac.params.beta_u','dnf_rch.params.h_u','dnf_rch.params.c_exc','dnf_rch.params.c_inh','dnf_rch.params.g_inh','dnf_rch.params.q_u',...
    'dnf_rch.params.beta_u'};
controlPosX = [0, 0, 1, 1, 0, 0, 1, 1, 2, 2, 2, 2] * controlFieldWidth;
controlPosY = [4, 3, 4, 3, 2, 1, 2, 1, 4, 3, 2, 1] * controlFieldHeight;
controlMin = [-10, 0, 0, 0, 0, 0, -10, 0, 0, 0, 0, 0];
controlMax = [0, 100, 100, 5, 1.5, 5.0, 0, 100, 100, 5, 1.5, 5.0];
textFormat = {'%0.1f', '%0.1f', '%0.1f', '%0.2f', '%0.2f', '%0.2f','%0.1f', '%0.1f', '%0.1f', '%0.2f', '%0.2f', '%0.2f'};

if debug==0
    close all
end

%%%%%%%%%%%%%%%%%%%%
% model parameters %
%%%%%%%%%%%%%%%%%%%%

dnf_SinputParams=initDNFParams();  %Spatial input field
dnf_SinputParams.c_inh = 30;
dnf_SinputParams.c_exc = 15;

dnf_sacParams=initDNFParams();    %Saccades field
dnf_sacParams.c_inh = 30;
dnf_sacParams.c_exc = 15;

dnf_rchParams=initDNFParams();    %Reach field
dnf_rchParams.c_inh = 30;
dnf_rchParams.c_exc = 15;

NcuedNeurons = 150;
Wcued        = normrnd(1,0.1,NcuedNeurons,dnf_rchParams.fieldSize);  %Weight matrix of context cue

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Stochastic optimal control framework
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
lqgParams_saccade=initLQGParamsSaccade();
lqgParams_reach  =initLQGParamsReach();

%%%%%%%%%%%%%%%%%%
% Protocols
%%%%%%%%%%%%%%%%%%
protocolParams=initSingleTargetProtocolParams();

%%%%%%%%%%%%%%%%%%%%%%%%%%
% simulation time course %
%%%%%%%%%%%%%%%%%%%%%%%%%%
simParams=initSimulationParams(protocolParams);


% set times at which field activities are stored (different variants)
% tStoreFields = [100, 200]; % select specific time steps
tStoreFields = 1:simParams.tMax; % store field activities at every time step

%%%%%%%%%%%%%%%%%%
% initialization %
%%%%%%%%%%%%%%%%%%
dnf_Sinput = initDNF(dnf_SinputParams, simParams.nTrials, tStoreFields);
dnf_sac    = initDNF(dnf_sacParams   , simParams.nTrials, tStoreFields);
dnf_rch    = initDNF(dnf_rchParams   , simParams.nTrials, tStoreFields);


stimulus_Sfield     = zeros(1,dnf_Sinput.params.fieldSize);
stimulus_sac        = zeros(1,dnf_sac.params.fieldSize);
stimulus_rch        = zeros(1,dnf_rch.params.fieldSize);
cue_context         = normrnd(0,0.05,1,NcuedNeurons);
effort_sac          = normrnd(0,0.05,1,dnf_sac.params.fieldSize);
effort_rch          = normrnd(0,0.05,1,dnf_rch.params.fieldSize);
CostSac =0;
CostRch =0;


%%%%%%%%%%%%%%%%%%%%
% Init Controllers %
%%%%%%%%%%%%%%%%%%%%
protocolTxt=uicontrol(fig,'Style','Text','Units','Norm','Position',[gapWidth, 5*controlFieldHeight, controlFieldWidth/2, controlFieldHeight], 'String', 'Protocol', 'HorizontalAlignment', 'left', 'BackgroundColor', 'w');
protocolDropdown=uicontrol(fig,'Style','Popupmenu', 'Units', 'Norm', 'Position', [gapWidth+.05, 5*controlFieldHeight, controlFieldWidth/2, controlFieldHeight], 'String', 'single target|two targets', 'Callback', @setProtocol);
cueTxt=uicontrol(fig,'Style','Text','Units','Norm','Position',[controlFieldWidth/2+2*gapWidth+.05, 5*controlFieldHeight, controlFieldWidth/2, controlFieldHeight], 'String', 'Cue', 'HorizontalAlignment', 'left', 'BackgroundColor', 'w');
cueDropdown=uicontrol(fig, 'Style', 'Popupmenu', 'Units', 'Norm', 'Position', [.075+controlFieldWidth/2+3*gapWidth, 5*controlFieldHeight, controlFieldWidth/2, controlFieldHeight], 'String', 'saccade|reach|free choice', 'Callback', @setCue);

startButton=uicontrol(fig, 'Style', 'Pushbutton', 'String', 'Start', 'Units', 'Norm', 'Position', [.075+controlFieldWidth+4*gapWidth, 5*controlFieldHeight+.0125, controlFieldWidth/4, .025], 'Callback', @startSim);
pauseButton=uicontrol(fig, 'Style', 'Pushbutton', 'String', 'Pause', 'Units', 'Norm', 'Position', [.075+controlFieldWidth+5*gapWidth+controlFieldWidth/4, 5*controlFieldHeight+.0125, controlFieldWidth/4, .025], 'Callback', @pauseSim);
stepButton=uicontrol(fig, 'Style', 'Pushbutton', 'String', 'Step', 'Units', 'Norm', 'Position', [.075+controlFieldWidth+6*gapWidth+2*controlFieldWidth/4, 5*controlFieldHeight+.0125, controlFieldWidth/4, .025], 'Callback', @stepSim);
stopButton=uicontrol(fig, 'Style', 'Pushbutton', 'String', 'Stop', 'Units', 'Norm', 'Position', [.075+controlFieldWidth+7*gapWidth+3*controlFieldWidth/4, 5*controlFieldHeight+.0125, controlFieldWidth/4, .025], 'Callback', @stopSim);

recordCheckbox=uicontrol(fig, 'Style', 'Checkbox', 'String', 'Record', 'Value', 0, 'Units', 'Norm', 'Position', [.075+controlFieldWidth+8*gapWidth+4*controlFieldWidth/4, 5*controlFieldHeight+.0125, controlFieldWidth/4, .025], 'Callback', @toggleRecord);
recordFilefield=uicontrol(fig, 'Style', 'Edit', 'String', 'DNF_LQG_simulation.avi', 'Units', 'Norm', 'Position', [.075+controlFieldWidth+9*gapWidth+5*controlFieldWidth/4, 5*controlFieldHeight+.0125, controlFieldWidth/2, .025]);

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


%%%%%%%%%%%%%%
% simulation %
%%%%%%%%%%%%%%
flag_evaluate_sac    = 1;
flag_evaluate_rch    = 1;
cnt_evaluate_sac     = 0;
cnt_evaluate_rch     = 0;
recordSim = 0;

% plot graphs
actPlot_u1=0;
outPlot_u1=0;
actPlot_u2=0;
outPlot_u2=0;
inPlot1=0;
inPlot2=0;
kernelPlot1=0;
kernelPlot2=0;
if debug>0
    
    axes(SpatialInAxes);
    cla;
    hold on;
    plot([0,dnf_Sinput.params.fieldSize-1],[0,0],'Linestyle',':','Linewidth',1);
    actPlot_u0 = plot(0:dnf_Sinput.params.fieldSize-1,dnf_Sinput.field_u,'color','b','Linewidth',3);
    outPlot_u0 = plot(0:dnf_Sinput.params.fieldSize-1,10*dnf_Sinput.output_u,'color','r','Linewidth',1);
    inPlot0 = plot(0:dnf_Sinput.params.fieldSize-1,stimulus_sac+dnf_Sinput.params.h_u,'color','g','Linewidth',1);
    set(gca,'ylim',[-15,15],'xlim',[0,dnf_Sinput.params.fieldSize-1],'Ytick',[-10,0,10]);
    ylabel('Spatial Input Field','Fontsize',12);
    hold off;
    
    
    axes(SacAxes);
    cla;
    hold on;
    plot([0,dnf_sac.params.fieldSize-1],[0,0],'Linestyle',':','Linewidth',1);
    actPlot_u1 = plot(0:dnf_sac.params.fieldSize-1,dnf_sac.field_u,'color','b','Linewidth',3);
    outPlot_u1 = plot(0:dnf_sac.params.fieldSize-1,10*dnf_sac.output_u,'color','r','Linewidth',1);
    inPlot1 = plot(0:dnf_sac.params.fieldSize-1,stimulus_sac+dnf_sac.params.h_u,'color','g','Linewidth',1);
    set(gca,'ylim',[-15,15],'xlim',[0,dnf_sac.params.fieldSize-1],'Ytick',[-10,0,10]);
    ylabel('saccade field','Fontsize',12);
    hold off;
    
    axes(ReachAxes);
    cla;
    hold on;
    plot([0,dnf_rch.params.fieldSize-1],[0,0],'Linestyle',':','Linewidth',1);
    actPlot_u2 = plot(0:dnf_rch.params.fieldSize-1,dnf_rch.field_u,'color','b','Linewidth',3);
    outPlot_u2 = plot(0:dnf_rch.params.fieldSize-1,10*dnf_rch.output_u,'color','r','Linewidth',1);
    inPlot2 = plot(0:dnf_rch.params.fieldSize-1,stimulus_rch+dnf_rch.params.h_u,'color','g','Linewidth',1);
    set(gca,'ylim',[-15,15],'xlim',[0,dnf_rch.params.fieldSize-1],'Ytick',[-10,0,10]);
    ylabel('reach field','Fontsize',12);
    hold off;
    
    axes(CueAxes);
    cla;
    hold on;
    plot([0,NcuedNeurons-1],[0,0],'Linestyle',':','Linewidth',1);
    outPlot_u3 = plot(0:NcuedNeurons-1,cue_context,'color','r','Linewidth',1);
    set(gca,'ylim',[-.1,1.1],'xlim',[0,NcuedNeurons-1],'Ytick',[-1,0,1]);
    ylabel('cue field','Fontsize',12);
    hold off;
    
    
    axes(kernelAxes);
    cla;
    hold on;
    plot([-dnf_sac.params.halfField,dnf_sac.params.halfField],[0,0],'Linestyle',':','Linewidth',1);
    kernelPlot1 = plot(-dnf_sac.params.halfField:dnf_sac.params.halfField, ...
        [zeros(1, dnf_sac.params.halfField-dnf_sac.kSize_uu) dnf_sac.kernel_uu zeros(1, dnf_sac.params.halfField-dnf_sac.kSize_uu)] - dnf_sac.params.g_inh, 'Color', 'r', 'Linewidth', 3);
    kernelPlot2 = plot(-dnf_rch.params.halfField:dnf_rch.params.halfField, ...
        [zeros(1, dnf_rch.params.halfField-dnf_rch.kSize_uu) dnf_rch.kernel_uu zeros(1, dnf_rch.params.halfField-dnf_rch.kSize_uu)] - dnf_rch.params.g_inh, 'Color', 'r', 'Linewidth', 3, 'LineStyle', '--');
    set(gca,'ylim',[-10,10],'xlim',[-dnf_rch.params.halfField,dnf_rch.params.halfField],'Ytick',[-10,-5,0,5,10]);
    ylabel('interaction kernel','Fontsize',12);
    hold off;
    
    axes(EffortSacAxes);
    cla;
    hold on;
    plot([0,dnf_sac.params.fieldSize-1],[0,0],'Linestyle',':','Linewidth',1);
    outPlot_u4 = plot(0:dnf_sac.params.fieldSize-1,effort_sac,'color','r','Linewidth',1);
    set(gca,'ylim',[-.2,1.5],'xlim',[0,dnf_sac.params.fieldSize-1],'Ytick',[-1,0,1]);
    ylabel('Effort Saccade','Fontsize',12);
    hold off;
    
    axes(EffortRchAxes);
    cla;
    hold on;
    plot([0,dnf_rch.params.fieldSize-1],[0,0],'Linestyle',':','Linewidth',1);
    outPlot_u5 = plot(0:dnf_rch.params.fieldSize-1,effort_sac,'color','r','Linewidth',1);
    set(gca,'ylim',[-.2,1.5],'xlim',[0,dnf_rch.params.fieldSize-1],'Ytick',[-1,0,1]);
    ylabel('Effort Reach','Fontsize',12);
    hold off;
    
end

cnt_movie = 0;

trial=1;

while true
	if simParams.running>0 || simParams.step>0
		% loop over trials
		while trial<=simParams.nTrials
	            while simParams.running<1 && simParams.step<1
                        pause(.001);
                    end
                    disp(['trial ' num2str(trial)]);
	    
                    protocolParams=resetProtocolParams(protocolParams);
                    simParams=resetSimParams(simParams, protocolParams);
                    flag_evaluate_sac=1;
                    flag_evaluate_rch=1;
                    cnt_evaluate_sac= 0;
                    cnt_evaluate_rch= 0;
		    dnf_Sinput = resetDNF(dnf_Sinput);
		    dnf_sac    = resetDNF(dnf_sac);
		    dnf_rch    = resetDNF(dnf_rch);
                    clear SaccadicMovement;
                    clear ReachingMovement;
                    
                    % loop over time steps
		    for t = 1 : simParams.tMax
		    	while simParams.running<1 && simParams.step<1
                            pause(.001);
                        end
                        if rem(t,20) == 0
	                    [MotorOutSac CostSac ActivePopSac] = LQGdnf(dnf_sac,simParams.output_u_sac_thr,...
                                lqgParams_saccade,simParams.Xeye,simParams.minEyeDist, ...
                                protocolParams.target_positions(1,:),protocolParams.target_positions(1,:),0);
	                    [MotorOutRch CostRch ActivePopRch] = LQGdnf(dnf_rch,simParams.output_u_rch_thr,...
                                lqgParams_reach,simParams.Xhand,simParams.minHandDist, ...
                                protocolParams.target_positions(1,:),protocolParams.target_positions(1,:),1);
	                end
			
                        if max(dnf_sac.output_u) > simParams.thr_perf_saccade
	                    if flag_evaluate_sac == 1 && simParams.flag_contact == 0 % Re-evaluate the saccadic policy
		                flag_evaluate_sac = 0;
		                cnt_evaluate_sac  = 1;
		                [SaccadicStateTotal, SaccadicMovement(:,t), CostSac] = ActionSelection(dnf_sac,...
                                    simParams.output_u_sac_thr,lqgParams_saccade,simParams.Xeye,...
                                    simParams.minEyeDist,protocolParams.target_positions(1,:),...
                                    protocolParams.target_positions(1,:),0);
		            else
		                cnt_evaluate_sac = cnt_evaluate_sac + 1;
		                SaccadicMovement(:,t) = SaccadicStateTotal([1,3],cnt_evaluate_sac);
		                simParams.Xeye(1:end-2) = SaccadicStateTotal(:,cnt_evaluate_sac);
                                simParams.eye_dist_targets=[];
                                for i=1:size(protocolParams.target_positions,1)
                                    rx=protocolParams.target_positions(i,1);
                                    ry=protocolParams.target_positions(i,2);
                                    simParams.eye_dist_targets=[simParams.eye_dist_targets; sqrt((SaccadicMovement(1,t)-rx)^2  + (SaccadicMovement(2,t)-ry)^2)];
                                end
                                simParams.minEyeDist=min(simParams.eye_dist_targets);
		                if cnt_evaluate_sac == simParams.revaluate_threshold
		                    flag_evaluate_sac = 1;
		                end
		            end
		            axes(BehaviorAxes);
		            saccPlot=plot(SaccadicMovement(1,t),SaccadicMovement(2,t),'o','color',[0.5 0.5 0.5]);
	                end
	
                        if max(dnf_rch.output_u) > simParams.thr_perf_saccade
	                    if flag_evaluate_rch == 1 && simParams.flag_contact == 0 % Re-evaluate the reaching policy
		                flag_evaluate_rch = 0;
		                cnt_evaluate_rch = 1;
		                [ReachingStateTotal, ReachingMovement(:,t), CostRch] = ActionSelection(dnf_rch,...
                                    simParams.output_u_rch_thr,lqgParams_reach,simParams.Xhand,...
                                    simParams.minHandDist,protocolParams.target_positions(1,:),...
                                    protocolParams.target_positions(1,:),1);
		            else
		                cnt_evaluate_rch = cnt_evaluate_rch + 1;
		                ReachingMovement(:,t) = ReachingStateTotal([1,2],cnt_evaluate_rch);
		                simParams.Xhand(1:end-2) = ReachingStateTotal(:,cnt_evaluate_rch);
                                simParams.hand_dist_targets=[];
                                for i=1:size(protocolParams.target_positions,1)
                                    rx=protocolParams.target_positions(i,1);
                                    ry=protocolParams.target_positions(i,2);
                                    simParams.hand_dist_targets=[simParams.hand_dist_targets; sqrt((ReachingMovement(1,t)-rx)^2  + (ReachingMovement(2,t)-ry)^2)];
                                end
                                simParams.minHandDist=min(simParams.hand_dist_targets);
		                if cnt_evaluate_rch == simParams.revaluate_threshold
		                    flag_evaluate_rch      = 1;
		                end
		            end
		            axes(BehaviorAxes);
		            reachPlot=plot(ReachingMovement(1,t),ReachingMovement(2,t),'ko');
	                end

                        if simParams.minEyeDist < simParams.Threshold_contact || simParams.minHandDist < simParams.Threshold_contact
                            simParams.flag_contact = 1;        
                        end

                        effort_sac      = normrnd(0,0.05,1,dnf_sac.params.fieldSize) + 10*CostSac;
                        effort_rch      = normrnd(0,0.05,1,dnf_rch.params.fieldSize) + 10*CostRch;
 
                        if protocolParams.protocol == 0 % single target is only presented			    
			    [protocolParams stimulus_Sfield stimulus_sac stimulus_rch cue_context]=runSingleTargetProtocol(t, simParams, protocolParams, lqgParams_saccade, lqgParams_reach, BehaviorAxes, dnf_Sinput, dnf_sac, dnf_rch, Wcued, effort_sac, effort_rch);
                        elseif protocolParams.protocol == 1 % two targets are only presented
                            [protocolParams stimulus_Sfield stimulus_sac stimulus_rch cue_context]=runTwoTargetProtocol(t, simParams, protocolParams, lqgParams_saccade, lqgParams_reach, BehaviorAxes, dnf_Sinput, dnf_sac, dnf_rch, Wcued, effort_sac, effort_rch);
                        end

                        disp(['t=' num2str(t)]);                        
			dnf_Sinput=runDNF(dnf_Sinput, stimulus_Sfield, tStoreFields, t);
			dnf_sac=runDNF(dnf_sac, stimulus_sac, tStoreFields, t);
			dnf_rch=runDNF(dnf_rch, stimulus_rch, tStoreFields, t);                        
			
			if debug>0
			    set(inPlot0, 'Ydata', stimulus_Sfield+dnf_Sinput.params.h_u);
			    set(actPlot_u0,'Ydata',dnf_Sinput.field_u);
			    set(outPlot_u0,'Ydata',10*dnf_Sinput.output_u);
			    set(inPlot1, 'Ydata', stimulus_sac+dnf_sac.params.h_u);
			    set(actPlot_u1,'Ydata',dnf_sac.field_u);
			    set(outPlot_u1,'Ydata',10*dnf_sac.output_u);
			    set(inPlot2, 'Ydata', stimulus_rch+dnf_rch.params.h_u);
			    set(actPlot_u2,'Ydata',dnf_rch.field_u);
			    set(outPlot_u2,'Ydata',10*dnf_rch.output_u);
			    set(outPlot_u3,'Ydata',cue_context);
			    set(outPlot_u4,'Ydata',effort_sac);
			    set(outPlot_u5,'Ydata',effort_rch);
			    drawnow;
			end
	
                        if simParams.step>0
                            simParams.step=0;
                        end
                        if recordSim>0
			    ff(t) = getframe(fig);
                        end
	
		    end
                    trial=trial+1;
		end
                if recordSim>0 && exist('ff')
                    disp('Writing movie file....');
                    recordFile=get(recordFilefield,'String');
		    movie2avi(ff, recordFile, 'compression', 'None');
                    disp('done!');
                    clear ff;
                end
	end
        pause(0.001);
end
return

nStoredFields = simParams.nTrials * length(tStoreFields);

figure;
subplot(3, 1, 1);
imagesc(dnf_sac.history_s');
xlabel('time');
ylabel('stimulus');
colorbar();
subplot(3, 1, 2);
imagesc(dnf_sac.history_u');
xlabel('time');
ylabel('activity u');
colorbar();
subplot(3, 1, 3);
imagesc(dnf_sac.history_output');
xlabel('time');
ylabel('output u');
colorbar();

figure;
subplot(3, 1, 1);
imagesc(dnf_rch.history_s');
xlabel('time');
ylabel('stimulus');
colorbar();
subplot(3, 1, 2);
imagesc(dnf_rch.history_u');
xlabel('time');
ylabel('activity u');
colorbar();
subplot(3, 1, 3);
imagesc(dnf_rch.history_output');
xlabel('time');
ylabel('output u');
colorbar();

% view evolution of field activities in each trial as mesh plot
nFieldsPerTrial = length(tStoreFields);
if 0
    disp('Press any key to iterate through trials');
    figure;
    for i = 1 : nStoredFields
        plot(0:dnf_sac.params.fieldSize-1, zeros(1, dnf_sac.params.fieldSize), ':k', ...
            1:dnf_sac.params.fieldSize, dnf_sac.history_s(i, :), '--g', ...
            1:dnf_sac.params.fieldSize, dnf_sac.history_mu(i, :) + dnf_sac.params.h_u, '-c', ...
            1:dnf_sac.params.fieldSize, dnf_sac.history_u(i, :), '-b');
        set(gca, 'XLim', [0 dnf_sac.params.fieldSize-1], 'YLim', [-15 15]);
        ylabel('activity u');
        drawnow;
        pause(0.01);
    end
    for i = 1 : nTrials
        subplot(2, 1, 1);
        mesh(1:dnf_sac.params.fieldSize, tStoreFields, ...
            dnf_sac.history_u((i-1)*nFieldsPerTrial+1 : i*nFieldsPerTrial, :));
        zlabel('activity u');
        subplot(2, 1, 2);
        mesh(1:dnf_sac.params.fieldSize, tStoreFields, ...
            dnf_sac.history_mu((i-1)*nFieldsPerTrial+1 : i*nFieldsPerTrial, :));
        zlabel('memory trace u');
        pause
    end
end

% view mesh plot of all stored field activities together
if 0
    figure;
    subplot(2, 1, 1);
    mesh(1:dnf_sac.params.fieldSize, 1:nStoredFields, dnf_sac.history_u(:, :));
    zlabel('activity u');
    subplot(2, 1, 2);
    mesh(1:dnf_sac.params.fieldSize, 1:nStoredFields, dnf_sac.history_mu(:, :));
    zlabel('memory trace u');
end

%update paramter values after slider changed
    function sliderCallback(hObject, eventdata) %#ok<INUSD>
        sliderChanged = find(hObject == sliders);
        
        paramName = controlParamNames{sliderChanged};
        tmp = get(sliders(sliderChanged), 'Value');
        set(textFields(sliderChanged), 'String', [paramName '=' num2str(tmp, textFormat{sliderChanged})]);
        eval([paramName '= tmp;']);
        dnf_sac.kernel_uu = dnf_sac.params.c_exc * gaussNorm(-dnf_sac.params.halfField:dnf_sac.params.halfField, 0, dnf_sac.params.sigma_exc) ...
            - dnf_sac.params.c_inh * gaussNorm(-dnf_sac.params.halfField:dnf_sac.params.halfField, 0, dnf_sac.params.sigma_inh) - dnf_sac.params.g_inh;
        dnf_sac.kernel_mu = dnf_sac.params.c_mu * gaussNorm(-dnf_sac.params.halfField:dnf_sac.params.halfField, 0, dnf_sac.params.sigma_mu);
        set(kernelPlot1,'YData',[zeros(1, dnf_sac.params.halfField-dnf_sac.kSize_uu) dnf_sac.kernel_uu zeros(1, dnf_sac.params.halfField-dnf_sac.kSize_uu)] - dnf_sac.params.g_inh);
        dnf_rch.kernel_uu = dnf_rch.params.c_exc * gaussNorm(-dnf_rch.params.halfField:dnf_rch.params.halfField, 0, dnf_rch.params.sigma_exc) ...
            - dnf_rch.params.c_inh * gaussNorm(-dnf_rch.params.halfField:dnf_rch.params.halfField, 0, dnf_rch.params.sigma_inh) - dnf_rch.params.g_inh;
        dnf_rch.kernel_mu = dnf_rch.params.c_mu * gaussNorm(-dnf_rch.params.halfField:dnf_rch.params.halfField, 0, dnf_rch.params.sigma_mu);
        set(kernelPlot2,'YData',[zeros(1, dnf_rch.params.halfField-dnf_rch.kSize_uu) dnf_rch.kernel_uu zeros(1, dnf_rch.params.halfField-dnf_rch.kSize_uu)] - dnf_rch.params.g_inh);
    end

    function setProtocol(hObject, eventdata)
        protocol=get(hObject, 'Value')-1;
        if protocol==0
            protocolParams=initSingleTargetProtocolParams();
        elseif protocol==1
            protocolParams=initTwoTargetProtocolParams();
        end
        simParams=initSimulationParams(protocolParams);
    end

    function setCue(hObject, eventdata)
        cueName=get(hObject, 'Value');
        protocolParams.cue=cueName-1;
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
        dnf_Sinput = resetDNF(dnf_Sinput);
        dnf_sac    = resetDNF(dnf_sac);
        dnf_rch    = resetDNF(dnf_rch);
        clear SaccadicMovement;
        clear ReachingMovement;
        trial=1;
        t=1;
        simParams.running=0;
    end

    function stepSim(hObject, eventdata)
        simParams.step=1;
    end

    function toggleRecord(hObject, eventdata)
        recordValue=get(hObject,'Value');
        recordSim=recordValue;
    end

% if max(dnf_sac.output_u)/max(dnf_rch.output_u) > 2
%     u_above_thr_id         = find(dnf_sac.output_u>output_u_sac_thr);
%     u_above_thr            = dnf_sac.output_u(u_above_thr_id)/sum(dnf_sac.output_u(u_above_thr_id));
%     NumOfActiveControllers = length(u_above_thr_id);
%     U_avg                  = 0;  %weighted average of the motor commands
%     current_distance       = sqrt((X1(1)-rx1)^2  + (X1(2)-ry1)^2); %current distance from target
%     Run_distance(t)        = current_distance ;%min(current_distance_T1,current_distance_T2);
%     disp(['Distance_to_contact=' num2str(Run_distance(t)]);
%     if Run_distance(t) < Threshold_contact
%         flag_contact = 1;
%         Contact_time = t;
%         tStimulusEnd = Contact_time; %When we contact the target, the stimuli disappear
%     end
%     for myControl = 1:NumOfActiveControllers
%         [x y] = pol2cart(u_above_thr_id(myControl)*pi/180, max(Run_distance_t1(t)));
%         x=x+X1(1); y=y+X1(2);
%         X1(end-1:end) = [x y];
%         [K,L,Cost,Xa,XSim,Xhat,CostSim,Unoise,Ufree] = kalman_lqg(A,B,C,C0, H,D,D0, E0, Q,R, X1,S1, NSim,Init,Niter );
%         [Cost_commands(myControl) Cost_accuracy(myControl) Overall_action_cost(myControl)] = ActionCost(Unoise,Q,R,Xhat,[rx1 ry1],[rx2 ry2]);
%         U_avg   = U_avg +  u_above_thr(myControl)*squeeze(Unoise);
%     end


end



