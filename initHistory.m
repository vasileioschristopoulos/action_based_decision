function history=initHistory(trials, freq, tMax, dnf_input_eye, dnf_input_hand, dnf_sac, dnf_rch, NCuedNeurons, NexpRewNeurons)

history.record_trials=[1:freq:trials];

n_recorded_trials=length(history.record_trials);
history.reaching_movements=zeros(n_recorded_trials,tMax,2);
history.saccade_movements=zeros(n_recorded_trials,tMax,2);

history.dnf_stim_eye=zeros(n_recorded_trials,tMax,dnf_input_eye.params.fieldSize);
history.dnf_stim_hand=zeros(n_recorded_trials,tMax,dnf_input_hand.params.fieldSize);
history.dnf_sac=zeros(n_recorded_trials,tMax,dnf_sac.params.fieldSize);
history.dnf_rch=zeros(n_recorded_trials,tMax,dnf_rch.params.fieldSize);
history.effort_rch=zeros(n_recorded_trials,tMax,dnf_rch.params.fieldSize);
history.effort_sac=zeros(n_recorded_trials,tMax,dnf_sac.params.fieldSize);
history.cue=zeros(n_recorded_trials,tMax,NCuedNeurons);

history.dopamine=zeros(1,trials);
history.cue_sac_w=zeros(trials,2,dnf_rch.params.fieldSize);
history.cue_rch_w=zeros(trials,2,dnf_rch.params.fieldSize);
history.expected_rew_w=zeros(trials,NexpRewNeurons,NexpRewNeurons);
