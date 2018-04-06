% Action cost for a SINGLE trajectory
function [Cost_commands Cost_accuracy Overall_action_cost] = ActionCost(u,R,X,TargetsLoc,effector,PrefDirect)

% Cost function: J = Xtend'*Q(:,:,tend)*Xtend + sum_{j=1}^N u(:,j)'*R*U(:,j)
%Inputs: u --> Sequence of control from current state to the goal
%        Q --> Q matrix of the cost function 
%        R --> R matrix of the cost function 
%        X --> State matrix (current state to the goal)
%        Target1 --> x,y coordinates of target 1
%        Target2 --> x,y coordinates of target 2

%Outputs: Cost_commands       --> Cost associated with the motor commands
%         Cost_accuracy       --> Cost associated with the accuracy at the end of the movement
%         Overall_action_cost --> Motor commands cost + accuracy cost
       
Cost_commands = 0;
TimeSteps = length(u); % Time step to reach the goal starting from the current state

for jj = 1:TimeSteps
    Cost_commands = Cost_commands + u(:,1,jj)'*R*u(:,1,jj);
end

X = squeeze(X);

Xeffector = X(1,end);
if effector == 1     % hand
    Yeffector = X(2,end);
else                 % eye
    Yeffector = X(3,end);
end

for myTarget = 1:size(TargetsLoc,1)
    direction_to_target(myTarget) = atan2d(TargetsLoc(myTarget,2)-Yeffector,TargetsLoc(myTarget,1)-Xeffector);
    distancet_to_target(myTarget) = sqrt((TargetsLoc(myTarget,1) - Xeffector)^2 + (TargetsLoc(myTarget,2) - Yeffector)^2);
    dev(myTarget)                 = abs(direction_to_target(myTarget) - PrefDirect);
end

Weight_cost = exp(-dev)./sum(exp(-dev));

Cost_accuracy = 0;
for myTarget   = 1:size(TargetsLoc,1)
    Cost_accuracy = Cost_accuracy + Weight_cost(myTarget)*distancet_to_target(myTarget);
end


Overall_action_cost = Cost_commands + 5*Cost_accuracy;

    
    
    
