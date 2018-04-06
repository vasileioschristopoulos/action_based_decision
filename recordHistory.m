function history=recordHistory(history, trial, reaching_movement, saccade_movement, dnf_input_eye, dnf_input_hand, dnf_sac, dnf_rch, cue, effort_sac, effort_rch, dopamine, Wcued, expected_reward)

if length(find(history.record_trials==trial))
    idx=find(history.record_trials==trial);
    history.reaching_movements(idx,:,:)=reaching_movement;
    history.saccade_movements(idx,:,:)=saccade_movement;
    history.dnf_stim_eye(idx,:,:)=dnf_input_eye.history_u;
    history.dnf_stim_hand(idx,:,:)=dnf_input_hand.history_u;
    history.dnf_sac(idx,:,:)=dnf_sac.history_u;
    history.dnf_rch(idx,:,:)=dnf_rch.history_u;
    history.cue(idx,:,:)=cue;
    history.effort_rch(idx,:,:)=effort_rch;
    history.effort_sac(idx,:,:)=effort_sac;
end
history.dopamine(trial)=dopamine;
history.cue_sac_w(trial,1,:)=mean(squeeze(Wcued(1,1:50,:)));
history.cue_sac_w(trial,2,:)=mean(squeeze(Wcued(1,51:100,:)));
history.cue_rch_w(trial,1,:)=mean(squeeze(Wcued(2,1:50,:)));
history.cue_rch_w(trial,2,:)=mean(squeeze(Wcued(2,51:100,:)));
history.expected_rew_w(trial,:,:)=expected_reward;
