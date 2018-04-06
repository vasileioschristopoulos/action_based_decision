function protocolParams=initCueFirstProtocolParams()

protocolParams=initProtocolParams();
protocolParams.protocol=1;

% set times of stimulus presentation
protocolParams.tStimulusStart = 50;
protocolParams.tTarget        = 100;
protocolParams.tStimulusEnd   = 350;

