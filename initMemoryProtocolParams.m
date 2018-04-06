function protocolParams=initMemoryProtocolParams()

protocolParams=initProtocolParams();
protocolParams.protocol=2;

% set times of stimulus presentation
protocolParams.tCueStart = 50;
protocolParams.tCueEnd   = 55;
protocolParams.tGo = 150;

protocolParams.check_go = 0;

