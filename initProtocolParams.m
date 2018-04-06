function protocolParams=initProtocolParams()

% Single target
protocolParams.protocol=0;
% 0: Instructed Saccade, 1:Instructed Reach, 2:Free Choice
protocolParams.cue=0;
% Whether to randomly show reach and saccade
protocolParams.random_cue=0;
protocolParams.stimCoeff=7.0;

protocolParams.check_draw_stimulus   = 0;
protocolParams.check_draw_context    = 0;
protocolParams.check_delete_stimulus = 0;

protocolParams.target_positions=[];

% 0=none, 1=left eye, 2=right eye, 3=left reach, 4=right reach
protocolParams.lesion=0;
