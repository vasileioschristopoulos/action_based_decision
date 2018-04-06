% simulate noisy trajectories
function [XSim] = SimTraj(A,B,C0,H,Q,Uavg,NSim,X1)

%Output: XSim, simulated state vector (i.e, position, velocity) without the
%target position (8x1)


% determine sizes
szX = size(A(1:end-2,1:end-2),1);
szY = size(H(:,1:end-2),1);
szC0 = size(C0,2);
N    = size(Q,3);

% Initialization of XSim
XSim = X1(1:4);


% if C0 ia scalar, set them to 0 matrices and adjust size
if length(C0(:))==1 & C0(1)==0,
    C0 = zeros(szX,1);
end


for k = 1:N-1
    Un          = Uavg(:,k);
    XSim(:,k+1) = A(1:end-2,1:end-2)*XSim(:,k) + B(1:end-2,:)*Un + C0*randn(szC0,NSim);
end
    


