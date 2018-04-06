function [protocolParams stim_input_eye stim_input_hand cue_context]=runCueFirstProtocol(t, simParams, protocolParams, dnf_input_eye, dnf_input_hand)
    
stim_input_eye=zeros(1,dnf_input_eye.params.fieldSize);
stim_input_hand=zeros(1,dnf_input_hand.params.fieldSize);
for i=1:size(protocolParams.target_positions,1)
    rx=protocolParams.target_positions(i,1);
    ry=protocolParams.target_positions(i,2);
    [theta,rho] = cart2pol(rx-simParams.xe,ry-simParams.ye);
    theta=(theta/pi)*180.0;
    stim_input_eye = stim_input_eye+protocolParams.stimCoeff*gauss(1:dnf_input_eye.params.fieldSize, round(theta), 5);
    [theta,rho] = cart2pol(rx-simParams.xh,ry-simParams.yh);
    theta=(theta/pi)*180.0;
    stim_input_hand = stim_input_hand+protocolParams.stimCoeff*gauss(1:dnf_input_hand.params.fieldSize, round(theta), 5);
end

cue_context   = [max(0,normrnd(0.05,0.05,1,50)), max(0,normrnd(0.05,0.05,1,50))];

if t<protocolParams.tStimulusStart
    stim_input_eye=zeros(1,dnf_input_eye.params.fieldSize);
    stim_input_hand=zeros(1,dnf_input_hand.params.fieldSize);
elseif t>=protocolParams.tStimulusStart && t<protocolParams.tTarget
    if protocolParams.cue == 0  %% SACCADE
        cue_context   = [normrnd(1,0.05,1,50), max(0,normrnd(0.05,0.05,1,50))];
        if protocolParams.check_draw_context == 0
            protocolParams.check_draw_context = 1;
            protocolParams.ContextCue1 = plot(0,21,'ro','MarkerSize',10,'MarkerFaceColor','r');
            set(gca,'ylim',[0,50],'xlim',[-45,45]);
            daspect([1 1 1])
            disp('cue saccade');
        end

    elseif protocolParams.cue == 1 %% REACHING
        cue_context   = [max(0,normrnd(0.05,0.05,1,50)), normrnd(1,0.05,1,50)];
        if protocolParams.check_draw_context == 0
            protocolParams.check_draw_context = 1;
            protocolParams.ContextCue2 = plot(0,19,'go','MarkerSize',10,'MarkerFaceColor','g');
            set(gca,'ylim',[0,50],'xlim',[-45,45]);
            daspect([1 1 1])
            disp('cue reaching');
        end

    else %% FREE CHOICE
        cue_context   = [max(0,normrnd(0.05,0.05,1,50)), max(0,normrnd(0.05,0.05,1,50))];
    end
    stim_input_eye=zeros(1,dnf_input_eye.params.fieldSize);
    stim_input_hand=zeros(1,dnf_input_hand.params.fieldSize);

elseif t>=protocolParams.tTarget
    if protocolParams.cue == 0  %% SACCADE
        cue_context   = [normrnd(1,0.05,1,50), max(0,normrnd(0.05,0.05,1,50))];
    elseif  protocolParams.cue == 1 %% REACHING
         cue_context   = [max(0,normrnd(0.05,0.05,1,50)), normrnd(1,0.05,1,50)];
    else %% FREE CHOICE
        cue_context   = [max(0,normrnd(0.05,0.05,1,50)), max(0,normrnd(0.05,0.05,1,50))];
    end        
    if protocolParams.check_draw_stimulus == 0
        protocolParams.check_draw_stimulus = 1;
        protocolParams.TargetLoc1=[];
        protocolParams.TargetLoc2=[];
        for i=1:size(protocolParams.target_positions,1)
            rx=protocolParams.target_positions(i,1);
            ry=protocolParams.target_positions(i,2);
            protocolParams.TargetLoc1 = [protocolParams.TargetLoc1; plot(rx,ry,'ko','MarkerSize',10,'MarkerFaceColor','k')];
            %protocolParams.TargetLoc2 = [protocolParams.TargetLoc2; plot(rx,ry,'ko','MarkerSize',70)];
            protocolParams.TargetLoc2 = [protocolParams.TargetLoc2; myCircle(simParams.Threshold_contact ,rx,ry,100,'k')];
        end
        set(gca,'ylim',[0,50],'xlim',[-45,45]);
        daspect([1 1 1])
    end
    if t>=protocolParams.tStimulusEnd
        stim_input_eye=zeros(1,dnf_input_eye.params.fieldSize);
        stim_input_hand=zeros(1,dnf_input_hand.params.fieldSize);
    end    
end

if simParams.flag_contact>0
    if protocolParams.check_delete_stimulus == 0
        for i=1:size(protocolParams.target_positions,1)
            delete(protocolParams.TargetLoc1(i))
            delete(protocolParams.TargetLoc2(i))
        end
        protocolParams.check_delete_stimulus = 1;
        protocolParams.tStimulusEnd = t;   %When we contact the target, the stimuli disappear
    end
end

