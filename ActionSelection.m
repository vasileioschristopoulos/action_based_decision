
function [EffectorStateTotal ControlX ControlY] = ActionSelection(MotorOut,ActivePop,dnf,lqgParams,Xstate,effector)

Weights    = dnf.output_u(ActivePop)/sum(dnf.output_u(ActivePop));
ControlX = 0;
ControlY = 0;
for myControls = 1:size(ActivePop,2)
    ControlX  = ControlX + Weights(myControls)*MotorOut(ActivePop(myControls)).x;
    ControlY  = ControlY + Weights(myControls)*MotorOut(ActivePop(myControls)).y;
end
EffectorStateTotal    = RunSimTraj(lqgParams,[ControlX ControlY]',Xstate,effector);

    
    


