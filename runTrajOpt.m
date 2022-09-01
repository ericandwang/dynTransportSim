close all
clear all
folder = '../dynTransportSim';
addpath(genpath(folder))

%% Options
% iteration loops
numIterations = 2;
% TOPP optimization
useObjectiveGradient = true;
useConstraintGradient = true;
useLinearization = false;
% animation
showAnimation = true;
showSnapshots = true;
movingWindow = false;

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
mu = 0.1/2;%0.1; %0.3 %0.5
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
s0 = [0; ...        % x_G -5
      0; ...        % dx_G
      0; ...        % y_G 5
      0; ...        % dy_G
      -pi/4; ...        % th
      0; ...        % dth
      r_GC; ...     % r_GC
      0; ...        % dr_GC
      0; ...
      0; ...        % th_C
      0];           % dth_C

% desired end position
x_des = [5; ... % x
         4; ... % y
         0];   % th
dx_des = [0; ... % dx
          0; ... % dy
          0];    % dth

%% Generate Friction Cone (object frame)
maxNormalForce = 10; % max normal force [N]
[fCone, vec] = generatefCone(param,maxNormalForce);

%% Spline Generation
controlPoints = 31;
porder = 5;
sBounds = [0 1];

cc = linspace(sBounds(1),sBounds(2),controlPoints);
xx = linspace(s0(1),x_des(1),controlPoints);
yy = linspace(s0(3),x_des(2),controlPoints);
%xx = (cc.^3).*x_des(1);
%yy = (2*cc - cc.^3)*x_des(2);

xp = spapi(optknt(cc,porder), cc, xx);
dxp = fnder(xp,1);
ddxp = fnder(xp,2);
yp = spapi(optknt(cc,porder), cc, yy);
dyp = fnder(yp,1);
ddyp = fnder(yp,2);

%% Iteration (TOPP <-> path)
times = zeros(numIterations,1);
for iii = 1:numIterations
%% TOPP Preoptimization
% equally spaced evaluation points and initial conditions
evalPoints = 31;
initialVel = 1;%1e-2;

if (iii == 1)
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
else
    P0 = psolve;
    dss0 = psolve(1:evalPoints);
    ddss0 = psolve(evalPoints+1:2*evalPoints);
    th0 = psolve(2*evalPoints+1:3*evalPoints);
    dth0 = psolve(3*evalPoints+1:4*evalPoints);
    ddth0 = psolve(4*evalPoints+1:5*evalPoints);
end

if (useConstraintGradient)
    disp('Calculating constraint gradient functions...')
    % calculate inequality constraint jacobian
    if (useLinearization)
        dcFun = dcGenTOPP_lin(r_GC, param, fCone, vec, dxp, ddxp, dyp, ddyp, ss0, evalPoints,th0);
    else
        dcFun = dcGenTOPP(r_GC, param, fCone, vec, dxp, ddxp, dyp, ddyp, ss0, evalPoints);
    end
    % calculate equality constraint jacobian
    dceqFun = dceqGenTOPP(ss0, evalPoints);
else
    dcFun = [];
    dceqFun = [];
end

%% TOPP Optimization

% objective function
fun = @(P) objTOPP(P,evalPoints);

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
Aeq = zeros(4,length(P0));
Aeq(1,1) = fnval(dxp,sBounds(1)); % velocity endpoints
Aeq(2,evalPoints) = fnval(dxp,sBounds(2));
Aeq(3,1) = fnval(dyp,sBounds(1));
Aeq(4,evalPoints) = fnval(dyp,sBounds(2));
Aeq(5,2*evalPoints+1) = 1; % angle endpoints
Aeq(6,3*evalPoints) = 1;
Aeq(7,4*evalPoints+1) = 1; % angular velocity endpoints
Aeq(8,5*evalPoints) = 1;
beq = [s0(2); dx_des(1); ...
       s0(4); dx_des(2); ...
       s0(5); x_des(3); ...
       s0(6); dx_des(3)];

% lower and upper bounds
lb = [zeros(evalPoints,1); ones(evalPoints,1).*-Inf; ...
      ones(evalPoints,1).*-Inf; ones(evalPoints,1).*-Inf; ...
      ones(evalPoints,1).*-Inf];
ub = [ones(evalPoints,1).*Inf; ones(evalPoints,1).*Inf; ...
      ones(evalPoints,1).*Inf; ones(evalPoints,1).*Inf; ...
      ones(evalPoints,1).*Inf];

% nonlinear constraint function
if (useLinearization)
    nonlcon = @(P) nonlinconTOPP_lin(P, s0, ss0, param, fCone, vec, tol, dxp, ddxp, dyp, ddyp, dcFun, dceqFun, th0);
else
    nonlcon = @(P) nonlinconTOPP(P, s0, ss0, param, fCone, vec, tol, dxp, ddxp, dyp, ddyp, dcFun, dceqFun);
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
    'SpecifyObjectiveGradient',useObjectiveGradient, ...
    'SpecifyConstraintGradient',useConstraintGradient, ...
    'MaxFunctionEvaluations', 5000, ...
    'Display','iter');

% invoke solver
psolve = fmincon(problem);
times(iii) = sum(1./(psolve(2:evalPoints-1)))*ss0(2)-ss0(1);

%% TOPP Results
P = reshape(psolve,numel(psolve)/5,5);
ss = ss0;
dss = P(:,1);
ddss = P(:,2);
ths = P(:,3);
dths = P(:,4);
ddths = P(:,5);

figure, plot(P(:,3))
title ('Orientation')
ylabel('angle [rad]')
figure, plot(P(:,1)), hold on, plot(P(:,2))
legend({'$\dot{s}$','$\ddot{s}$'},'Interpreter','latex')
title('$\dot{s}$ and $\ddot{s}$','Interpreter','latex')
%figure, plot(fnval(xp,ss0), fnval(yp,ss0))

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

%% Animation
if (showAnimation)
    % initialize animation
    CoR_radius = 0.005;
    figure(1000); clf;
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
            [objectH/2, objectH/2, -objectH/2, -objectH/2, objectH/2]+y_G(i)+r_GC(2), 'color', [0.8500    0.3250    0.0980]);
        CoR = rectangle('Position',[x_G(i)-CoR_radius/2 y_G(i)-CoR_radius/2 ...
            CoR_radius CoR_radius],'Curvature',[1 1]);
        rotate(hand,[0 0 1],rad2deg(th(i)),[x_G(i),y_G(i),0])
        rotate(object,[0 0 1],rad2deg(th(i)),[x_G(i),y_G(i),0])
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
else
    plot3(p(1,:),p(2,:),p(3,:),'r')
end
%hold on
%scatter3(pSave(1,1),pSave(2,1),pSave(3,1),'k')
%plot3(pSave(1,:),pSave(2,:),pSave(3,:),'b')
%legend({'','','','','','\gamma = 0.1','','\gamma = 0.4'})
xlabel('$\ddot{x}_{c} [m/s^2]$','interpreter','latex')
ylabel('$\ddot{y}_{c} [m/s^2]$','interpreter','latex')
zlabel('$\ddot{\theta}_{c} [rad/s^2]$','interpreter','latex')
title('Gravito-Inertial Wrench Constraints')

%% Repathing

endPoints = [xx(1) xx(end) yy(1) yy(end)];
xxyy0 = [xx(2:end-1) yy(2:end-1)];

% objective function
fun = @(xxyy) objPATH(P, r_GC, param, fCone, vec, ss0, cc, endPoints, xxyy, porder, controlPoints);

% lower and upper bounds
splinePoints = controlPoints-2;
xxmin = min(endPoints(1),endPoints(2));
xxmax = max(endPoints(1),endPoints(2));
yymin = min(endPoints(3),endPoints(4));
yymax = max(endPoints(3),endPoints(4));
lb = [ones(splinePoints,1).*xxmin; ones(splinePoints,1).*yymin];
ub = [ones(splinePoints,1).*xxmax; ones(splinePoints,1).*yymax];

% nonlinear constraint function
nonlcon = @(xxyy) nonlinconPATH(P, r_GC, param, fCone, vec, ss0, cc, endPoints, xxyy, porder, controlPoints, tol);

% optimization
problem.x0 = xxyy0;
problem.objective = fun;
problem.Aineq = [];
problem.bineq = [];
problem.Aeq = [];
problem.beq = [];
problem.lb = lb;
problem.ub = ub;
problem.nonlcon = nonlcon;
problem.solver = 'fmincon';
problem.options = optimoptions('fmincon', 'Display','iter');

% invoke solver
disp('Replanning path...')
xxyysolve = fmincon(problem);

%% Path results
% save old path
xxold = xx;
yyold = yy;
xpold = xp;
ypold = yp;
dxpold = dxp;
dypold = dyp;
ddxpold = ddxp;
ddypold = ddyp;

% new path
xx = [endPoints(1) xxyysolve(1:splinePoints) endPoints(2)];
yy = [endPoints(3) xxyysolve(splinePoints+1:end) endPoints(4)];
xp = spapi(optknt(cc,porder), cc, xx);
dxp = fnder(xp,1);
ddxp = fnder(xp,2);
yp = spapi(optknt(cc,porder), cc, yy);
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
else
    plot3(p(1,:),p(2,:),p(3,:),'r')
end
%hold on
%scatter3(pSave(1,1),pSave(2,1),pSave(3,1),'k')
%plot3(pSave(1,:),pSave(2,:),pSave(3,:),'b')
%legend({'','','','','','\gamma = 0.1','','\gamma = 0.4'})
xlabel('$\ddot{x}_{c} [m/s^2]$','interpreter','latex')
ylabel('$\ddot{y}_{c} [m/s^2]$','interpreter','latex')
zlabel('$\ddot{\theta}_{c} [rad/s^2]$','interpreter','latex')
title('Gravito-Inertial Wrench Constraints')

% Show difference in path
figure, plot(fnval(xpold,ss0), fnval(ypold,ss0))
hold on, plot(fnval(xp,ss0), fnval(yp,ss0))
end