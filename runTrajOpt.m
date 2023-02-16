close all
clear all
folder = '../dynTransportSim';
addpath(genpath(folder))

%% Options
% iteration loops
numIterations = 15;
% convergence rate for transit time
alpha = 0.2; % 0.35 for example good
beta = 0; % thought is to constrain how much replanning can push points to the edge of the gravito-inertial wrench boundary
% Gen files
genFiles = 0;
% Warm start?
warmStart = 1;
% TOPP optimization
useTOPPObjectiveGradient = true;
useTOPPConstraintGradient = true;
useLinearization = false;
% PATH optimization
usePATHObjectiveGradient = false;
usePATHConstraintGradient = false;
convexConeApproximation = true;
% animation
showAnimation = true;
showSnapshots = 6;
movingWindow = false;
showForces = true;
xForceStretch = 3;

%% Parameter Setup
% hand and object dimensions
handW = 0.26;
handH = 0.01;
objectW = 0.10;
objectH = 0.10*1;
% parameters
m_G = 1.2; % kg
I_G = 0.1; % kg*m^2
m_C = 0.4; % kg
I_C = 0.05; % kg*m^2
r_p1C = [objectW/2; objectH/2]; % m
r_p2C = [-objectW/2; objectH/2]; % m
g = 9.81;
mu = 0.1/2*5;%0.1; %0.3 %0.5
param = [m_G; I_G; m_C; I_C; r_p1C; r_p2C; g; mu; handW; handH; objectW; objectH];
% controller actuation limits (bidirectional)
u_max = [30; 30; 20];
% tolerance for distance from plane constraints
tol = 0.05;%0.01 (0.05 = conservative);
% Parameters for control barrier function cone
a = (2.4/25)^2; % most conservative (2.5/25)^2
b = (8.9/25)^2; % most conservative (9/25)^2

% initial conditions
r_GC = [0.07; objectH/2+handH/2];
s0 = [-5; ...        % x_G -5
      0; ...        % dx_G
      5; ...        % y_G 5
      0; ...        % dy_G
      -3/4*pi; ...        % th (0)
      0; ...        % dth
      r_GC; ...     % r_GC
      0; ...        % dr_GC
      0; ...
      0; ...        % th_C
      0];           % dth_C
% s0 = [-5; ...        % x_G -5
%       0; ...        % dx_G
%       0; ...        % y_G 5
%       0; ...        % dy_G
%       0; ...        % th (0)
%       0; ...        % dth
%       r_GC; ...     % r_GC
%       0; ...        % dr_GC
%       0; ...
%       0; ...        % th_C
%       0];           % dth_C

% desired end position
x_des = [5; ... % x
         1; ... % y
         -pi/4];   % th
dx_des = [1; ... % dx
          1; ... % dy
          0];    % dth
% x_des = [0; ... % x
%          0; ... % y
%          0];   % th
% dx_des = [0; ... % dx
%           0; ... % dy
%           0];    % dth
s_des = [x_des(1); dx_des(1); x_des(2); dx_des(2); x_des(3); dx_des(3); ...
         s0(7:end)];

%% Generate Friction Cone (object frame)
maxNormalForce = 10; % max normal force [N]
[fCone, vec] = generatefCone(param,maxNormalForce);

%% Spline Generation
collPoints = 15;
porder = 3+1;
sBounds = [0 1];

cc = linspace(sBounds(1),sBounds(2),collPoints);
xx = linspace(s0(1),x_des(1),collPoints);
yy = linspace(s0(3),x_des(2),collPoints);

knotVec = [ones(1,porder-1).*sBounds(1) linspace(sBounds(1),sBounds(2),collPoints) ones(1,porder-1).*sBounds(2)];
interpPoints = linspace(sBounds(1),sBounds(2),collPoints+porder-2);
xp = spapi(knotVec, interpPoints, interp1(cc,xx,interpPoints));
dxp = fnder(xp,1);
ddxp = fnder(xp,2);
yp = spapi(knotVec, interpPoints, interp1(cc,yy,interpPoints));
%yp = spmak(yp.knots,[0 ones(1, length(yp.coefs)-2) 0]); % CCC remove immediately
%dyp = fnder(yp,1);
%ddyp = fnder(yp,2);
[xp, yp] = intermediatePlanBridge(s0(1), s0(3), x_des(1), x_des(2), knotVec, porder); % CCC remove immediately
dxp = fnder(xp,1);
ddxp = fnder(xp,2);
yp.coefs(porder:end-porder+1) = 1;
dyp = fnder(yp,1);
ddyp = fnder(yp,2);


%% Dynamically Feasible Spline Preinitialization
if (warmStart)
    disp('Calculating dynamically feasible warm start...')
    numPoints = 100;
    % statically infeasible IC to statically feasible point (0)
    [th_0, dth_0, ddth_0, tTotal_0, xp_0, yp_0] = intermediatePlanDynamic(s0, x_des, fCone, param, knotVec, porder, 1);
    % statically feasible point to statically infeasible FC (4)
    [th_4, dth_4, ddth_4, tTotal_4, xp_4, yp_4] = intermediatePlanDynamic(s_des, [s0(1); s0(3); s0(5)], fCone, param, knotVec, porder, -1);
    dxp_4 = fnder(xp_4,1); ddxp_0 = fnder(xp_4,2);
    dyp_4 = fnder(yp_4,1); ddyp_0 = fnder(yp_4,2);
    
    % Static path (1)
    [xp_1, yp_1] = intermediatePlanStatic(xp_0, yp_0, knotVec, porder, 1);
    
    % Static path (3)
    [xp_3, yp_3] = intermediatePlanStatic(xp_4, yp_4, knotVec, porder, -1);
    
    % Static path (2)
    %x2 = linspace(fnval(xp_1,1),fnval(xp_3,0),numPoints);
    %y2 = linspace(fnval(yp_1,1),fnval(yp_3,0),numPoints);
    %s2 = linspace(0,1,length(x2));
    %xp_2 = spap2(knotVec, porder, s2, x2);
    %yp_2 = spap2(knotVec, porder, s2, y2);
    [xp_2, yp_2] = intermediatePlanBridge(fnval(xp_1,1), fnval(yp_1,1), fnval(xp_3,0), fnval(yp_3,0), knotVec, porder);
    
    % Offsetting path parameter s
    xp_1.knots = xp_1.knots + 1;
    yp_1.knots = yp_1.knots + 1;
    xp_2.knots = xp_2.knots + 2;
    yp_2.knots = yp_2.knots + 2;
    xp_3.knots = xp_3.knots + 3;
    yp_3.knots = yp_3.knots + 3;
    xp_4.knots = xp_4.knots + 4;
    yp_4.knots = yp_4.knots + 4;
    
    % Stitching splines together
    xSplines = {xp_0, xp_1, xp_2, xp_3, xp_4};
    ySplines = {yp_0, yp_1, yp_2, yp_3, yp_4};
    numSplines = length(xSplines);
    sBounds = [sBounds(1) sBounds(2)*numSplines];
    collPoints = collPoints*numSplines; % CCC can increase to + 7 or + 5
    knotVec = [ones(1,porder-1).*sBounds(1) linspace(sBounds(1),sBounds(2),collPoints) ones(1,porder-1).*sBounds(2)];
    numPoints = length(knotVec) - (porder-1)*2;
    s = linspace(sBounds(1),sBounds(2),numPoints);
    ddx = zeros(1,length(s));
    ddy = zeros(1,length(s));
    for i = 1:numSplines
        ddx = ddx + fnval(fnder(xSplines{i},2),s);
        ddy = ddy + fnval(fnder(ySplines{i},2),s);
    end
    s = linspace(sBounds(1),sBounds(2),length(ddx));
    ddxp = spap2(knotVec(3:end-2), porder-2, s, ddx);
    ddyp = spap2(knotVec(3:end-2), porder-2, s, ddy);
    dxp = fnint(ddxp,fnval(fnder(xSplines{1},1),0));
    dyp = fnint(ddyp,fnval(fnder(ySplines{1},1),0));
    xp = fnint(dxp,fnval(xSplines{1},0));
    yp = fnint(dyp,fnval(ySplines{1},0));
    
    % debugging intermediate path termination CCC
    x_des = [fnval(xp,sBounds(2)); ... % x
             fnval(yp,sBounds(2)); ... % y
             0];   % th
    dx_des = [0; ... % dx
              0; ... % dy
              0];    % dth
end

%% TOPP/PATH Generating Analytical Gradient Functions
% number of constraint evaluation points
evalPoints = (collPoints-1)*3+1;
ss0 = linspace(sBounds(1),sBounds(2),evalPoints)';
% object frame y acceleration limit
accelLim = 100;

% pre-calculation for sparse gradient assembler
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

if (genFiles)
    disp('Calculating TOPP constraint gradient functions...')
    %dcGenTOPP3(r_GC, param, fCone, vec, evalPoints, accelLim);
    dceqGenTOPP2(ss0, evalPoints);

    disp('Calculating PATH constraint gradient functions...')
    xcoefs0 = xp.coefs';
    ycoefs0 = yp.coefs';
    coefs0 = [xcoefs0; ycoefs0];
    numCoefs = length(coefs0);
    %dcGenPATH2()
    disp('Calculating PATH objective gradient functions...')
    [basisVectors, constraintSlopes] = changeOfBasis(fCone);
    %objGenPATH_convex2()
end

%% Iteration (TOPP <-> path)
times = zeros(numIterations,1);
for iii = 1:numIterations
%% TOPP Preoptimization
% equally spaced evaluation points and initial conditions
initialVel = 1e-2; %1;

if (iii == 1)
    if (warmStart)
        % CCC replaced naive initialization with feasible spline
        % preinitialization results
        dss0 = ones(evalPoints,1).*1/tTotal_0;
        dss0(ss0 >= 4) = 1/tTotal_4;    
        ddss0 = zeros(evalPoints,1);
        for ii = 1:length(ddss0)-1
            ddss0(ii) = (dss0(ii+1)^2 - dss0(ii)^2)/(2*ss0(ii+1)-ss0(ii));
        end
        th0 = th_0(ss0*tTotal_0);
        th0(ss0 >= 4) = th_4((ss0(ss0>=4)-4)*tTotal_4);
        dth0 = dth_0(ss0*tTotal_0);
        dth0(ss0 >= 4) = dth_4((ss0(ss0>=4)-4)*tTotal_4);
        ddth0 = ddth_0(ss0*tTotal_0);
        ddth0(ss0 >= 4) = ddth_4((ss0(ss0>=4)-4)*tTotal_4);
        P0 = [dss0; ddss0; th0; dth0; ddth0];
    else
        ss0 = linspace(sBounds(1),sBounds(2),evalPoints)';
        dss0 = ones(evalPoints,1).*initialVel;
        ddss0 = zeros(evalPoints,1);
        for ii = 1:length(ddss0)-1
            ddss0(ii) = (dss0(ii+1)^2 - dss0(ii)^2)/(2*ss0(ii+1)-ss0(ii));
        end
        th0 = ones(evalPoints,1).*s0(5);
        dth0 = zeros(evalPoints,1);
        ddth0 = zeros(evalPoints,1);
        P0 = [dss0; ddss0; th0; dth0; ddth0];
    end
else
    P0 = psolve;
    dss0 = psolve(1:evalPoints);
    ddss0 = psolve(evalPoints+1:2*evalPoints);
    th0 = psolve(2*evalPoints+1:3*evalPoints);
    dth0 = psolve(3*evalPoints+1:4*evalPoints);
    ddth0 = psolve(4*evalPoints+1:5*evalPoints);
end

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
    dceqFun = @(states) dceqFunGen(states);
else
    dcFun = [];
    dceqFun = [];
end

%% TOPP Optimization

% objective function
fun = @(P) objTOPP(P,evalPoints);
if (iii > 1)
    tt_prev = fun(P0);
end

% linear inequality constraints
% monotonic path parameter (Note: removed but can reuse for variable delta
% path lengths)
Aineq1 = diag(ones(evalPoints,1),0);
Aineq2 = diag(-1.*ones(evalPoints,1),1);
Aineq = Aineq1(1:evalPoints-1,1:evalPoints) + ...
    Aineq2(1:evalPoints-1,1:evalPoints);
Aineq = [Aineq zeros(size(Aineq,1),evalPoints*5)];
bineq = zeros(evalPoints-1,1);

% linear equality constraints
Aeq = zeros(8,length(P0));
Aeq(1,1) = fnval(dxp,sBounds(1)); % velocity endpoints
Aeq(2,evalPoints) = fnval(dxp,sBounds(2));
Aeq(3,1) = fnval(dyp,sBounds(1));
Aeq(4,evalPoints) = fnval(dyp,sBounds(2));
Aeq(5,2*evalPoints+1) = 1; % angle endpoints
Aeq(6,3*evalPoints) = 1;
Aeq(7,3*evalPoints+1) = 1; % angular velocity endpoints
Aeq(8,4*evalPoints) = 1;
beq = [s0(2); dx_des(1); ...
       s0(4); dx_des(2); ...
       s0(5); x_des(3); ...
       s0(6); dx_des(3)];

% CCC removing terminal velocity constraints
%selectedRows = [1 3 5 6 7];
%Aeq = Aeq(selectedRows,:);
%beq = beq(selectedRows);

% lower and upper bounds
lb = [zeros(evalPoints,1); ones(evalPoints,1).*-Inf; ...
      ones(evalPoints,1).*-Inf; ones(evalPoints,1).*-Inf; ...
      ones(evalPoints,1).*-Inf];
ub = [ones(evalPoints,1).*Inf; ones(evalPoints,1).*Inf; ...
      ones(evalPoints,1).*Inf; ones(evalPoints,1).*Inf; ...
      ones(evalPoints,1).*Inf];

% nonlinear constraint function
if (useLinearization) % CCC not aligning with original constraint function now (decomissioned)
    nonlcon = @(P) nonlinconTOPP_lin(P, s0, ss0, param, fCone, vec, tol, dxp, ddxp, dyp, ddyp, dcFun, dceqFun, th0);
else
    if (iii == numIterations)
        nonlcon = @(P) nonlinconTOPP(P, s0, ss0, param, fCone, vec, tol, xp, dxp, ddxp, yp, dyp, ddyp, dcFun, dceqFun, accelLim);
    else
        %nonlcon = @(P) nonlinconTOPP(P, s0, ss0, param, fCone, vec, tol*(10-(iii-1)/2), xp, dxp, ddxp, yp, dyp, ddyp, dcFun, dceqFun, accelLim);
        if (iii >= 100) % 5 works CCC
            nonlcon = @(P) nonlinconTOPP(P, s0, ss0, param, fCone, vec, tol*5, xp, dxp, ddxp, yp, dyp, ddyp, dcFun, dceqFun, accelLim);
        else
            nonlcon = @(P) nonlinconTOPP(P, s0, ss0, param, fCone, vec, tol, xp, dxp, ddxp, yp, dyp, ddyp, dcFun, dceqFun, accelLim);
        end
    end
end

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
    'MaxFunctionEvaluations', 3000/30, ...
    'Display','iter');

% invoke solver
psolve = fmincon(problem);
tt_max = fun(psolve);
tts(iii) = tt_max;

%% Max lambda TOPP
% finding transit time constraint based on convergence rate
if (iii == 1)
    tt_constrain = alpha*tt_max + (1-alpha)*tt_max*2;
else
    tt_constrain = alpha*tt_max + (1-alpha)*tt_prev;
end

% augmenting state matrix with lambda
P0 = [psolve; 0.01];

% objective function
fun = @(P) objTOPPLambda(P);

% linear inequality constraints
% monotonic path parameter (Note: removed but can reuse for variable delta
% path lengths)
Aineq1 = diag(ones(evalPoints,1),0);
Aineq2 = diag(-1.*ones(evalPoints,1),1);
Aineq = Aineq1(1:evalPoints-1,1:evalPoints) + ...
    Aineq2(1:evalPoints-1,1:evalPoints);
Aineq = [Aineq zeros(size(Aineq,1),evalPoints*5+1)]; % making sure to include lambda
bineq = zeros(evalPoints-1,1);

% linear equality constraints
Aeq = zeros(8,length(P0));
Aeq(1,1) = fnval(dxp,sBounds(1)); % velocity endpoints
Aeq(2,evalPoints) = fnval(dxp,sBounds(2));
Aeq(3,1) = fnval(dyp,sBounds(1));
Aeq(4,evalPoints) = fnval(dyp,sBounds(2));
Aeq(5,2*evalPoints+1) = 1; % angle endpoints
Aeq(6,3*evalPoints) = 1;
Aeq(7,3*evalPoints+1) = 1; % angular velocity endpoints
Aeq(8,4*evalPoints) = 1;
beq = [s0(2); dx_des(1); ...
       s0(4); dx_des(2); ...
       s0(5); x_des(3); ...
       s0(6); dx_des(3)];

% CCC removing terminal velocity constraints
%selectedRows = [1 3 5 6 7];
%Aeq = Aeq(selectedRows,:);
%beq = beq(selectedRows);

% lower and upper bounds
lb = [zeros(evalPoints,1); ones(evalPoints,1).*-Inf; ...
      ones(evalPoints,1).*-Inf; ones(evalPoints,1).*-Inf; ...
      ones(evalPoints,1).*-Inf; 0];
ub = [ones(evalPoints,1).*Inf; ones(evalPoints,1).*Inf; ...
      ones(evalPoints,1).*Inf; ones(evalPoints,1).*Inf; ...
      ones(evalPoints,1).*Inf; Inf];

% nonlinear constraint function
nonlcon = @(P) nonlinconTOPPLambda(P, s0, ss0, param, fCone, vec, tol, xp, dxp, ddxp, yp, dyp, ddyp, dcFun, dceqFun, accelLim, tt_constrain);

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
    'MaxFunctionEvaluations', 3000/30, ...
    'Display','iter');

% invoke solver
psolve = fmincon(problem);
% removing lambda from state vector
lambda = psolve(end);
lambdas(iii) = lambda;
psolve = psolve(1:end-1);


%% Save Transit Time
%times(iii) = sum(1./(psolve(2:evalPoints-1)))*ss0(2)-ss0(1);
times(iii) = sum(2./(psolve(1:evalPoints-1,1)+psolve(2:evalPoints,1)))*ss0(2)-ss0(1);

%% TOPP Results
P = reshape(psolve,numel(psolve)/5,5);
ss = ss0;
dss = P(:,1);
ddss = P(:,2);
ths = P(:,3);
dths = P(:,4);
ddths = P(:,5);

if (iii == numIterations)
    figure, plot(P(:,3))
    title ('Orientation')
    ylabel('angle [rad]')
    figure, plot(P(:,1)), hold on, plot(P(:,2))
    legend({'$\dot{s}$','$\ddot{s}$'},'Interpreter','latex')
    title('$\dot{s}$ and $\ddot{s}$','Interpreter','latex')
end
%figure, plot(fnval(xp,ss0), fnval(yp,ss0))

%% Friction Cone

figure(100)
hold on
patch([0; fCone(1,1); fCone(1,2)],[0;fCone(2,1);fCone(2,2)],[0;fCone(3,1);fCone(3,2)],[0.5,0,0.5],'FaceAlpha',0.1)
patch([0; fCone(1,2); fCone(1,3)],[0;fCone(2,2);fCone(2,3)],[0;fCone(3,2);fCone(3,3)],[0.5,0,0.5],'FaceAlpha',0.1)
patch([0; fCone(1,3); fCone(1,4)],[0;fCone(2,3);fCone(2,4)],[0;fCone(3,3);fCone(3,4)],[0.5,0,0.5],'FaceAlpha',0.1)
patch([0; fCone(1,4); fCone(1,1)],[0;fCone(2,4);fCone(2,1)],[0;fCone(3,4);fCone(3,1)],[0.5,0,0.5],'FaceAlpha',0.1)
p = zeros(3,length(ss));
for i = 1:length(ss)
    ds = P(i,1);
    dds = P(i,2);
    th = P(i,3);
    dth = P(i,4);
    ddth = P(i,5);
    ax_G = fnval(ddxp,ss(i))*ds^2 + fnval(dxp,ss(i))*dds;
    ay_G = fnval(ddyp,ss(i))*ds^2 + fnval(dyp,ss(i))*dds;
    g = param(9);
    R_GC = norm(r_GC);
    th_GC = angle(r_GC(1) + 1i*r_GC(2));
    p(:,i) = [cos(th) sin(th) -R_GC*sin(th_GC); ...
              -sin(th) cos(th) R_GC*cos(th_GC); ...
              0 0 1]*[ax_G; ay_G; ddth] + ...
              [g*sin(th)-R_GC*dth^2*cos(th_GC); ...
              g*cos(th)-R_GC*dth^2*sin(th_GC); 0];
end
%for i = 1:length(t)
%    pSave(:,i) = worldToObject(t(i),sS(i,:)',param,uS(i,:));
%end
scatter3(p(1,1),p(2,1),p(3,1),'k');
if (iii == 1)
    plot3(p(1,:),p(2,:),p(3,:),'b')
elseif (iii == 2)
    plot3(p(1,:),p(2,:),p(3,:),'r')
else
    plot3(p(1,:),p(2,:),p(3,:),'g')
end
%hold on
%scatter3(pSave(1,1),pSave(2,1),pSave(3,1),'k')
%plot3(pSave(1,:),pSave(2,:),pSave(3,:),'b')
%legend({'','','','','','\gamma = 0.1','','\gamma = 0.4'})
xlabel('$\ddot{x}_{c} [m/s^2]$','interpreter','latex')
ylabel('$\ddot{y}_{c} [m/s^2]$','interpreter','latex')
zlabel('$\ddot{\theta}_{c} [rad/s^2]$','interpreter','latex')
title('Gravito-Inertial Wrench Constraints')

%% Animation
% forward time scaling
dT = 2.*(ss(2:end)-ss(1:end-1))./(dss(2:end)+dss(1:end-1));
tss = cumsum([0; dT]);
dt = min(dT)/10;
dx = fnval(dxp,ss).*dss;
dy = fnval(dyp,ss).*dss;
dth = dths;
t = 0:dt:tss(end);
dx = interp1(tss,dx,t);
dy = interp1(tss,dy,t);
dth = interp1(tss,dth,t);
x_G = cumsum(dx.*dt)' + s0(1);
y_G = cumsum(dy.*dt)' + s0(3);
th = cumsum(dth.*dt)' + s0(5);
forces = individualForces(vec,p);
forces = forces./max(forces).*0.9*objectH;
fLx = forces(1,:); fLy = forces(2,:);
fRx = forces(3,:); fRy = forces(4,:);
fLx = interp1(tss,fLx,t); fLy = interp1(tss,fLy,t);
fRx = interp1(tss,fRx,t); fRy = interp1(tss,fRy,t);
fLx = fLx*xForceStretch; fRx = fRx*xForceStretch;

if (showAnimation)
    % initialize animation
    CoR_radius = 0.005;
    figure(1002); clf;
    axis equal
    
    % show animation
    for i = 1:1+29*showSnapshots:length(t)
        if (showSnapshots)
            hold on
        else
            hold off
        end
        hand = plot([-handW/2, handW/2, handW/2, -handW/2, -handW/2]+x_G(i), ...
            [handH/2, handH/2, -handH/2, -handH/2, handH/2]+y_G(i), 'color', [0    0.4470    0.7410]);
        hold on
        object = plot([-objectW/2, objectW/2, objectW/2, -objectW/2, -objectW/2]+x_G(i)+r_GC(1), ...
            [objectH/2, objectH/2, -objectH/2, -objectH/2, objectH/2]+y_G(i)+r_GC(2), 'color', 'black'); % [0.8500    0.3250    0.0980]
        CoR = rectangle('Position',[x_G(i)-CoR_radius/2 y_G(i)-CoR_radius/2 ...
            CoR_radius CoR_radius],'Curvature',[1 1]);
        rotate(hand,[0 0 1],rad2deg(th(i)),[x_G(i),y_G(i),0])
        rotate(object,[0 0 1],rad2deg(th(i)),[x_G(i),y_G(i),0])
        if (showForces)
            forceL = plot([-objectW/2 -objectW/2+fLx(i)]+x_G(i)+r_GC(1), ...
                [-objectH/2 -objectH/2+fLy(i)]+y_G(i)+r_GC(2), 'color', 'red');
            forceR = plot([objectW/2 objectW/2+fRx(i)]+x_G(i)+r_GC(1), ...
                [-objectH/2 -objectH/2+fRy(i)]+y_G(i)+r_GC(2), 'color', 'red');
            frictionConeL = plot([-objectW/2 -objectW/2-objectH/2*mu*xForceStretch -objectW/2 -objectW/2+objectH/2*mu*xForceStretch]+x_G(i)+r_GC(1), ...
                [-objectH/2 0 -objectH/2 0]+y_G(i)+r_GC(2), '--', 'color', 'green');
            frictionConeR = plot([objectW/2 objectW/2-objectH/2*mu*xForceStretch objectW/2 objectW/2+objectH/2*mu*xForceStretch]+x_G(i)+r_GC(1), ...
                [-objectH/2 0 -objectH/2 0]+y_G(i)+r_GC(2), '--', 'color', 'green');
            rotate(forceL,[0 0 1],rad2deg(th(i)),[x_G(i),y_G(i),0])
            rotate(forceR,[0 0 1],rad2deg(th(i)),[x_G(i),y_G(i),0])
            rotate(frictionConeL,[0 0 1],rad2deg(th(i)),[x_G(i),y_G(i),0])
            rotate(frictionConeR,[0 0 1],rad2deg(th(i)),[x_G(i),y_G(i),0])
        end
        if (movingWindow)
            axis equal
            axis manual
            plot(x_G(1:i),y_G(1:i))
        else
            plot(x_G(1:i),y_G(1:i),'--')
            xlim([min([x_G;y_G])-1 max([x_G;y_G])+1])
            ylim([min([x_G;y_G])-1 max([x_G;y_G])+1])
        end
        ylabel('y [m]')
        xlabel('x [m]')
        drawnow
    end
end

%% PATH Preoptimization

% initial B-spline control points
xcoefs0 = xp.coefs';
ycoefs0 = yp.coefs';
coefs0 = [xcoefs0; ycoefs0];
numCoefs = length(coefs0);
[basis, dbasis, ddbasis, posEndpoints, velEndpoints] = splineBasisCoefs(knotVec, porder);

if (usePATHConstraintGradient)
    disp('Calculating PATH constraint gradient functions...')
    % calculate inequality constraint jacobian
    dcFun = dcGenPATH(P, r_GC, param, fCone, vec, ss0, knotVec, numCoefs, accelLim);
    % calculate equality constraint jacobian
    dceqFun = dceqGenPATH(numCoefs, velEndpoints, s0, dx_des);
else
    dcFun = [];
    dceqFun = [];
end

if (convexConeApproximation)
    [basisVectors, constraintSlopes] = changeOfBasis(fCone);
end

if (usePATHObjectiveGradient)
    disp('Calculating PATH objective gradient functions...')
    % calculate inequality constraint jacobian
    if (convexConeApproximation)
        gradfFun = objGenPATH_convex(P, r_GC, param, ss0, knotVec, numCoefs, basisVectors, constraintSlopes);
    else
        gradfFun = objGenPATH(P, r_GC, param, fCone, vec, ss0, knotVec, numCoefs);
    end
    % calculate equality constraint jacobian
else
    gradfFun = [];
end

%% Repathing PATH

% objective function
if (convexConeApproximation)
    fun = @(coefs) objPATH_convex(P, r_GC, param, ss0, knotVec, coefs, basisVectors, constraintSlopes, gradfFun);
else
    fun = @(coefs) objPATH(P, r_GC, param, fCone, vec, ss0, knotVec, coefs, gradfFun);
end

% control point bounds
%lb = ones(numCoefs,1).*-10;
%ub = ones(numCoefs,1).*10;

% endpoint constraints and acceleration constraints
Aeq = zeros(4,numCoefs);
beq = zeros(4,1);
% position endpoints
Aeq(1,1:numCoefs/2) = posEndpoints(1,:);
Aeq(2,1:numCoefs/2) = posEndpoints(2,:);
Aeq(3,numCoefs/2+1:end) = posEndpoints(1,:);
Aeq(4,numCoefs/2+1:end) = posEndpoints(2,:);
beq(1) = s0(1); beq(2) = x_des(1); beq(3) = s0(3); beq(4) = x_des(2);


% add in constraints for velocity endpoints IMPORTANT CCC
% accelerations CCC not implemented
%Aeq2 = [ddImpulses' zeros(size(ddImpulses')); ...
%        zeros(size(ddImpulses')) ddImpulses'; ...
%        -ddImpulses' zeros(size(ddImpulses')); ...
%        zeros(size(ddImpulses')) -ddImpulses'];
%numAccel = size(ddImpulses,2);
%limAccel = 1000;
%beq2 = [ones(numAccel*2,1).*limAccel; ones(numAccel*2,1).*-limAccel];
%Aeq = [Aeq1];
%beq = [beq1];
% object frame y acceleration limit
%accelLim = max(p(2,:))*1.2;


% nonlinear constraint function
nonlcon = @(coefs) nonlinconPATH(P, r_GC, param, fCone, vec, ss0, knotVec, coefs, tol + beta*lambda, accelLim, velEndpoints, s0, dx_des, dcFun, dceqFun);

% optimization
problem.x0 = coefs0;
problem.objective = fun;
problem.Aineq = [];
problem.bineq = [];
problem.Aeq = Aeq;
problem.beq = beq;
problem.lb = [];
problem.ub = [];
problem.nonlcon = nonlcon;
problem.solver = 'fmincon';
problem.options = optimoptions('fmincon', ...
    'SpecifyObjectiveGradient',usePATHObjectiveGradient, ...
    'SpecifyConstraintGradient',usePATHConstraintGradient, ... 
    'Display','iter');

% invoke solver
disp('Replanning path...')
coefs = fmincon(problem);

%% Path results
% save old path
xpold{iii} = xp;
ypold{iii} = yp;
dxpold{iii} = dxp;
dypold{iii} = dyp;
ddxpold{iii} = ddxp;
ddypold{iii} = ddyp;

% new path
numCoefs = length(coefs);
xcoefs = coefs(1:numCoefs/2);
ycoefs = coefs(numCoefs/2+1:end);
xp = spmak(knotVec,xcoefs');
dxp = fnder(xp,1);
ddxp = fnder(xp,2);
yp = spmak(knotVec,ycoefs');
dyp = fnder(yp,1);
ddyp = fnder(yp,2);

% new friction cone
figure(101)
hold on
patch([0; fCone(1,1); fCone(1,2)],[0;fCone(2,1);fCone(2,2)],[0;fCone(3,1);fCone(3,2)],[0.5,0,0.5],'FaceAlpha',0.1)
patch([0; fCone(1,2); fCone(1,3)],[0;fCone(2,2);fCone(2,3)],[0;fCone(3,2);fCone(3,3)],[0.5,0,0.5],'FaceAlpha',0.1)
patch([0; fCone(1,3); fCone(1,4)],[0;fCone(2,3);fCone(2,4)],[0;fCone(3,3);fCone(3,4)],[0.5,0,0.5],'FaceAlpha',0.1)
patch([0; fCone(1,4); fCone(1,1)],[0;fCone(2,4);fCone(2,1)],[0;fCone(3,4);fCone(3,1)],[0.5,0,0.5],'FaceAlpha',0.1)
p = zeros(3,length(ss));
for i = 1:length(ss)
    ds = P(i,1);
    dds = P(i,2);
    th = P(i,3);
    dth = P(i,4);
    ddth = P(i,5);
    ax_G = fnval(ddxp,ss(i))*ds^2 + fnval(dxp,ss(i))*dds;
    ay_G = fnval(ddyp,ss(i))*ds^2 + fnval(dyp,ss(i))*dds;
    g = param(9);
    R_GC = norm(r_GC);
    th_GC = angle(r_GC(1) + 1i*r_GC(2));
    p(:,i) = [cos(th) sin(th) -R_GC*sin(th_GC); ...
              -sin(th) cos(th) R_GC*cos(th_GC); ...
              0 0 1]*[ax_G; ay_G; ddth] + ...
              [g*sin(th)-R_GC*dth^2*cos(th_GC); ...
              g*cos(th)-R_GC*dth^2*sin(th_GC); 0];
end
%for i = 1:length(t)
%    pSave(:,i) = worldToObject(t(i),sS(i,:)',param,uS(i,:));
%end
scatter3(p(1,1),p(2,1),p(3,1),'k');
if (iii == 1)
    plot3(p(1,:),p(2,:),p(3,:),'b')
elseif (iii == 2)
    plot3(p(1,:),p(2,:),p(3,:),'r')
else
    plot3(p(1,:),p(2,:),p(3,:),'g')
end
%hold on
%scatter3(pSave(1,1),pSave(2,1),pSave(3,1),'k')
%plot3(pSave(1,:),pSave(2,:),pSave(3,:),'b')
%legend({'','','','','','\gamma = 0.1','','\gamma = 0.4'})
xlabel('$\ddot{x}_{c} [m/s^2]$','interpreter','latex')
ylabel('$\ddot{y}_{c} [m/s^2]$','interpreter','latex')
zlabel('$\ddot{\theta}_{c} [rad/s^2]$','interpreter','latex')
title('Gravito-Inertial Wrench Constraints')

end

%% Post Iteration Results
% Show difference in path
splot = linspace(ss0(1),ss0(end),length(ss0)*5);
figure
for i = 1:iii
    hold on
    plot(fnval(xpold{i},splot), fnval(ypold{i},splot))
end
plot(fnval(xp,splot), fnval(yp,splot))
xlabel('x [m]')
ylabel('y [m]')
title('Path Iteration')

% Transit time plot
figure, plot(times)
xlim([1 length(times)])
ylabel('transit time [sec]')
xlabel('iteration #')
