function [psolve, tt_prev] = runTOPP(s0, tol, evalPoints, r_GC, param, fCone, vec, accelLim, ss0, dceqFunFile, initialVel, ...
    xp, yp, useTOPPObjectiveGradient, useTOPPConstraintGradient, beq, sBounds)

    % spline
    dxp = fnder(xp,1);
    ddxp = fnder(xp,2);
    dyp = fnder(yp,1);
    ddyp = fnder(yp,2);

    % gradients
    dcFun1 = dcGenTOPP3(r_GC, param, fCone, vec, 1, accelLim);
    selectionMatrix_ = eye(evalPoints*5);
    selectionMatrixL = selectionMatrix_;
    for i = 1:5 % iterate over number of states
        selectionMatrixL(1+evalPoints*(i-1):evalPoints*i,:) = ...
            selectionMatrix_(i:5:end,:);
    end
    selectionMatrixR = selectionMatrix_;
    nCon = 5; % number of constraints
    selectionMatrixR(:,evalPoints*(nCon-1)+1:end) = selectionMatrix_(:,nCon:nCon:end);
    for i = 1:evalPoints
        selectionMatrixR(:,1+(nCon-1)*(i-1):(nCon-1)*i) = selectionMatrix_(:,1+nCon*(i-1):nCon*(i-1)+nCon-1);
    end
    dss0 = ones(evalPoints,1).*initialVel;
    ddss0 = zeros(evalPoints,1);
    for ii = 1:length(ddss0)-1
        ddss0(ii) = (dss0(ii+1)^2 - dss0(ii)^2)/(2*ss0(ii+1)-ss0(ii));
    end
    th0 = zeros(evalPoints,1);
    dth0 = zeros(evalPoints,1);
    ddth0 = zeros(evalPoints,1);
    P0 = [dss0; ddss0; th0; dth0; ddth0];
    if (useTOPPConstraintGradient)
        % calculate inequality constraint jacobian
        %dcFun = dcGenTOPP(r_GC, param, fCone, vec, xp, dxp, ddxp, yp, dyp, ddyp, ss0, evalPoints, accelLim);
        dxp_ = fnval(dxp,ss0);
        ddxp_ = fnval(ddxp,ss0);
        dyp_ = fnval(dyp,ss0);
        ddyp_ = fnval(ddyp,ss0);
        splines = [dxp_; ddxp_; dyp_; ddyp_];
        %dcFun = @(states) dcFunGen(states, splines);
        dcFun = @(states) dcGenTOPPAssemble(states, splines, dcFun1, selectionMatrixL, selectionMatrixR);
        % calculate equality constraint jacobian
        %dceqFun = dceqGenTOPP(ss0, evalPoints);
        dceqFun = @(states) dceqFunFile(states);
    else
        dcFun = [];
        dceqFun = [];
    end

    % objective function
    fun = @(P) objTOPP(P,evalPoints);
    tt_prev = fun(P0);

    % linear equality constraints
    Aeq = zeros(8,evalPoints*5);
    Aeq(1,1) = fnval(fnder(xp,1),sBounds(1)); % velocity endpoints
    Aeq(2,evalPoints) = fnval(fnder(xp,1),sBounds(2));
    Aeq(3,1) = fnval(fnder(yp,1),sBounds(1));
    Aeq(4,evalPoints) = fnval(fnder(yp,1),sBounds(2));
    Aeq(5,2*evalPoints+1) = 1; % angle endpoints
    Aeq(6,3*evalPoints) = 1;
    Aeq(7,3*evalPoints+1) = 1; % angular velocity endpoints
    Aeq(8,4*evalPoints) = 1;
    
    % lower and upper bounds
    lb = [zeros(evalPoints,1); ones(evalPoints,1).*-Inf; ...
          ones(evalPoints,1).*-Inf; ones(evalPoints,1).*-Inf; ...
          ones(evalPoints,1).*-Inf];
    ub = [ones(evalPoints,1).*Inf; ones(evalPoints,1).*Inf; ...
          ones(evalPoints,1).*Inf; ones(evalPoints,1).*Inf; ...
          ones(evalPoints,1).*Inf];
    
    % nonlinear constraint function
    nonlcon = @(P) nonlinconTOPP(P, s0, ss0, param, fCone, vec, tol, xp, dxp, ddxp, yp, dyp, ddyp, dcFun, dceqFun, accelLim);

    % optimization
    problem.x0 = P0;
    problem.objective = fun;
    problem.Aineq = [];
    problem.bineq = [];
    problem.Aeq = Aeq;
    problem.beq = beq;
    problem.lb = lb;
    problem.ub = ub;
    problem.nonlcon = nonlcon;
    problem.solver = 'fmincon';
    problem.options = optimoptions('fmincon', ...
        'SpecifyObjectiveGradient',useTOPPObjectiveGradient, ...
        'SpecifyConstraintGradient',useTOPPConstraintGradient, ...
        'Display','iter');
    % invoke solver
    psolve = fmincon(problem);
end

