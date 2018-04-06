function protocolParams=initTargetFirstProtocolParams()

protocolParams=initProtocolParams();
protocolParams.protocol=0;

% set times of stimulus presentation
protocolParams.tStimulusStart = 50;
protocolParams.tCue           = 100;
protocolParams.tStimulusEnd   = 350;

