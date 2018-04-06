% simulate noisy trajectories
function [XSim] = RunSimTraj(lqgParams,Uavg,X1,effector)

%Output: XSim, simulated state vector (i.e, position, velocity) without the
%target position --> hand: (8x1) and eye: 4x1) 


% determine sizes
szX  = size(lqgParams.A(1:end-2,1:end-2),1);
szY  = size(lqgParams.H(:,1:end-2),1);
szC0 = size(lqgParams.C0,2);
N    = size(lqgParams.Q,3);

% Initialization of XSim
if effector == 1     % hand
    XSim = X1(1:8);
else                % eye
    XSim = X1(1:4);
end


% if C0 ia scalar, set them to 0 matrices and adjust size
if length(lqgParams.C0(:))==1 & lqgParams.C0(1)==0,
    lqgParams.C0 = zeros(szX,1);
end


for k = 1:N-1
    Un          = Uavg(:,k);
    XSim(:,k+1) = lqgParams.A(1:end-2,1:end-2)*XSim(:,k) + lqgParams.B(1:end-2,:)*Un + lqgParams.C0*randn(szC0,lqgParams.NSim);
end
    


