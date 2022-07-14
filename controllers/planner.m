function [uSave, sSave, s0Save, tSave] = planner(controllerBandwidth, tDuration, s0, param, fCone, vec, tol, jacob)
% nonlinear trajectory optimization

% states
th = s0(5);

% parameters
g = param(9);

nHorizon = tDuration*controllerBandwidth;
dT = 1/controllerBandwidth;

costVector = [0; 0; 0; 0; 0; 0; 1; 1; 1];
%finalCostVector = [100; 1; 100; 1; 10; 1; 1; 1; 1];
fun = @(x) objgrad(x,nHorizon,costVector,costVector);
x0 = [s0(1:6,1);zeros(3+9*(nHorizon-1),1)];
nonlincon_ = @(x) nonlincon2grad(nHorizon,x,s0,param,fCone,vec,tol,jacob); % change for right nonlincon function
A = [];
b = [];
% linear equality constraints
Aeq = zeros(6*nHorizon + 6,9*nHorizon);
beq = zeros(6*nHorizon + 6,1);
% copy matrix to be pasted into Aeq
mCopy = [1 dT 0 0  0 0  1/2*dT^2 0 0 -1 0 0 0 0 0; ...
         0 1  0 0  0 0  dT 0 0       0 -1 0 0 0 0; ...
         0 0  1 dT 0 0  0 1/2*dT^2 0 0 0 -1 0 0 0; ...
         0 0  0 1  0 0  0 dT 0        0 0 0 -1 0 0; ...
         0 0  0 0  1 dT 0 0 1/2*dT^2 0 0 0 0 -1 0; ...
         0 0  0 0  0 1  0 0 dT       0 0 0 0 0  -1];
for i = 1:nHorizon+1
    if (i == 1)
        Aeq(1:6,1:6) = eye(6);
    elseif(i == nHorizon+1)
        Aeq(1+6*(i-1):6+6*(i-1), ...
            1+9*(i-2):6+9*(i-2)) = eye(6);
    else
        Aeq(1+6*(i-1):6+6*(i-1), ...
            1+9*(i-2):1+9*(i-2)+size(mCopy,2)-1) = mCopy;
    end
end
beq(1:6,1) = s0(1:6,1);
%bCopy = [30; 30; 30; 30; pi; 30; 30; 30; 30];
bCopy = [inf; inf; inf; inf; pi; inf; inf; inf; inf];
lb = repmat(-bCopy,[nHorizon,1]);
ub = repmat(bCopy,[nHorizon,1]);

options = optimoptions('fmincon','MaxFunctionEvaluations',5000,'SpecifyObjectiveGradient',true, ...
    'SpecifyConstraintGradient',true');
xSolve = fmincon(fun,x0,A,b,Aeq,beq,lb,ub,nonlincon_,options);


u = xSolve(7:9);

end

