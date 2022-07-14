function u = controller(dT, nHorizon, s, param, fCone, vec, tol, jacob, a_des, x_g)
% MPC controller that minimizes cost function over receding horizon

% states
th = s(5);

% parameters
g = param(9);

% Converting desired acceleration to object frame
%rot = @(th) [cos(th) -sin(th); sin(th) cos(th)];
%a_des_ = [rot(-th)*a_des(1:2); a_des(3)];

if (nHorizon == 1)
    fun = @(x) norm(x(7:9)-a_des);
    x0 = [s(1:6,1);a_des];
else
    costVector = [0; 0; 0; 0; 0; 0; 1; 1; 1];
    finalCostVector = [100; 1; 100; 1; 10; 1; 1; 1; 1];
    %costMatrix = diag([repmat(costVector,[nHorizon-1,1]); finalCostVector]);
    %fun = @(x) x'*costMatrix*x;
    fun = @(x) objgrad(x,x_g,nHorizon,costVector,finalCostVector);
    x0 = [s(1:6,1);zeros(3+9*(nHorizon-1),1)];
end
nonlincon_ = @(x) nonlincon2grad(nHorizon,x,s,param,fCone,vec,tol,jacob); % change for right nonlincon function
A = [];
b = [];
% linear equality constraints
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
%bCopy = [30; 30; 30; 30; pi; 30; 30; 30; 30];
bCopy = [inf; inf; inf; inf; pi; inf; inf; inf; inf];
lb = repmat(-bCopy,[nHorizon,1]);
ub = repmat(bCopy,[nHorizon,1]);

options = optimoptions('fmincon','MaxFunctionEvaluations',5000,'SpecifyObjectiveGradient',true, ...
    'SpecifyConstraintGradient',true');
if (nHorizon == 1)
    nonlincon_ = @(x) nonlincon2(nHorizon,x,s,param,fCone,vec,tol);
    xSolve = fmincon(fun,x0,A,b,Aeq,beq,lb,ub,nonlincon_);    
else
    xSolve = fmincon(fun,x0,A,b,Aeq,beq,lb,ub,nonlincon_,options);
end

u = xSolve(7:9);

end

