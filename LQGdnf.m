%function [MotorOut CostField u_above_thr_planning_id]  = LQGdnf(dnf,output_u_threshold,lqgParams,Xstate,distance,Target1,Target2,effector)
function [MotorOut CostField u_above_thr_planning_id]  = LQGdnf(dnf,output_u_threshold,lqgParams,Xstate,distance,TargetsLoc,effector)
    CostField   = zeros(1,dnf.params.fieldSize);  
    MotorOut(dnf.params.fieldSize).x  = [];
    MotorOut(dnf.params.fieldSize).y  = [];
    
    u_above_thr_planning_id = find(dnf.output_u>output_u_threshold);
    if isempty(u_above_thr_planning_id) == 0 %Planning
        for jj = 1:length(u_above_thr_planning_id)
            [x y] = pol2cart(u_above_thr_planning_id(jj)*pi/180, distance);
            if effector == 1     % hand
                x=x+Xstate(1); y=y+Xstate(2);
            else   % eye
                x=x+Xstate(1); y=y+Xstate(3);
            end 
            Xstate(end-1:end) = [x y];
            [K,L,Cost,Xa,XSim,Xhat,CostSim,Unoise,Ufree] = runLQG(lqgParams,Xstate);
            %[Cost_command(jj) Cost_accuracy(jj) Overall_action_cost(jj)] = ActionCost(Unoise,lqgParams.Q,lqgParams.R,Xhat,Target1,Target2);
            [Cost_command(jj) Cost_accuracy(jj) Overall_action_cost(jj)] = ActionCost(Unoise,lqgParams.R,Xhat,TargetsLoc,effector,u_above_thr_planning_id(jj));
            CostField(u_above_thr_planning_id(jj)) = Overall_action_cost(jj);
            MotorOut(u_above_thr_planning_id(jj)).x =  squeeze(Unoise(1,1,:));
            MotorOut(u_above_thr_planning_id(jj)).y =  squeeze(Unoise(2,1,:));
        end
        CostField = CostField./sum(CostField);
    end
end