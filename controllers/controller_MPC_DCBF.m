function u = controller_MPC_DCBF(dT, nHorizon, s, param, Bp, Bp_initial, gamma, x_g)
% MPC controller that minimizes cost function over receding horizon

% states
th = s(5);

% parameters
g = param(9);

costVector = [0; 0; 0; 0; 0; 1; 1; 1; 1];
finalCostVector = [100; 1; 100; 1; 10; 1; 1; 1; 1];
fun = @(x) objgrad(x,x_g,nHorizon,costVector,finalCostVector);
x0 = [s(1:6,1);zeros(3+9*(nHorizon-1),1)];

nonlincon_ = @(x) nonlinconPlaneCBF(nHorizon,x,s,param,Bp,Bp_initial,gamma);
A = [];
b = [];
% linear equality constraints (dynamics)
Aeq = zeros(6*nHorizon,9*nHorizon);
beq = zeros(6*nHorizon,1);
% copy matrix to be pasted into Aeq
mCopy = [1 dT 0 0  0 0  1/2*dT^2 0 0 -1 0 0 0 0 0; ...
         0 1  0 0  0 0  dT 0 0       0 -1 0 0 0 0; ...
         0 0  1 dT 0 0  0 1/2*dT^2 0 0 0 -1 0 0 0; ...
         0 0  0 1  0 0  0 dT 0        0 0 0 -1 0 0; ...
         0 0  0 0  1 dT 0 0 1/2*dT^2 0 0 0 0 -1 0; ...
         0 0  0 0  0 1  0 0 dT       0 0 0 0 0  -1];
for i = 1:nHorizon
    if (i == 1)
        Aeq(1:6,1:6) = eye(6);
    else
        Aeq(1+6*(i-1):6+6*(i-1), ...
            1+9*(i-2):1+9*(i-2)+size(mCopy,2)-1) = mCopy;
    end
end
beq(1:6,1) = s(1:6,1);
bCopy = [inf; inf; inf; inf; pi; inf; inf; inf; inf];
lb = repmat(-bCopy,[nHorizon,1]);
ub = repmat(bCopy,[nHorizon,1]);

options = optimoptions('fmincon','MaxFunctionEvaluations',5000,'SpecifyObjectiveGradient',true, ...
    'SpecifyConstraintGradient',false,'Algorithm','sqp');
%x0 = fmincon(fun,x0,A,b,Aeq,beq,lb,ub);
xSolve = fmincon(fun,x0,A,b,Aeq,beq,lb,ub,nonlincon_,options);

% Debug printing out CBF trajectories
% figure(1001)
% Bvec_ = zeros(4*(nHorizon+1),1);
% Bvec_(1:4,1) = Bp_initial;
% for i = 1:nHorizon
%     s_u = xSolve(1+9*(i-1):9+9*(i-1));
%     Bvec_(1+4*i:4+4*i,1) = Bp(s_u);
% end
% plot(Bvec_(1:4:end))
% hold on, plot(Bvec_(2:4:end)), plot(Bvec_(3:4:end)), plot(Bvec_(4:4:end)), hold off

u = xSolve(7:9);

end

