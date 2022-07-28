close all
clear all
folder = '../dynTransportSim';
addpath(genpath(folder))

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
mu = 0.7;%0.1; %0.3 %0.5
param = [m_G; I_G; m_C; I_C; r_p1C; r_p2C; g; mu; handW; handH; objectW; objectH];
% controller actuation limits (bidirectional)
u_max = [30; 30; 20];
% tolerance for distance from plane constraints
tol = 0.05;%0.01 (0.05 = conservative);
% Parameters for control barrier function cone
a = (2.4/25)^2; % most conservative (2.5/25)^2
b = (8.9/25)^2; % most conservative (9/25)^2

% initial conditions
s0 = [0; ... % x_G -5
      0; ...   % dx_G
      0; ...   % y_G 5
      0; ...   % dy_G
      0; ...   % th
      0; ...   % dth
      0.07; ...   % r_GC 0.07
      objectH/2+handH/2; ...
      0; ... %dr_GC
      0; ...
      0; ... %th_C
      0]; %dth_C

% desired end position
x_des = [5; ... % x
         4; ... % y
         0;];   % th
dx_des = [0; ... % dx
          0; ... % dy
          0];    % dth

%% Generate Friction Cone (object frame)
maxNormalForce = 10; % max normal force [N]
[fCone, vec] = generatefCone(param,maxNormalForce);

%% Spline Generation
controlPoints = 11;
porder = 5;
sBounds = [0 1];

cc = linspace(sBounds(1),sBounds(2),controlPoints);
xx = linspace(s0(1),x_des(1),controlPoints);
yy = linspace(s0(3),x_des(3),controlPoints);
%xx = (cc.^3).*x_des(1);
%yy = (2*cc - cc.^3)*x_des(2);

xp = spapi(optknt(cc,porder), cc, xx);
dxp = fnder(xp,1);
ddxp = fnder(xp,2);
yp = spapi(optknt(cc,porder), cc, yy);
dyp = fnder(yp,1);
ddyp = fnder(yp,2);

%% TOPP
% equally spaced evaluation points and initial conditions
evalPoints = 31;
ss0 = linspace(sBounds(1),sBounds(2),evalPoints)';
dss0 = ones(evalPoints,1);
ddss0 = zeros(evalPoints,1);
for ii = 1:length(ddss0)-1
    ddss0(ii) = (dss0(ii+1)^2 - dss0(ii)^2)/(2*ss0(ii+1)-ss0(ii));
end
th0 = zeros(evalPoints,1);
dth0 = zeros(evalPoints,1);
ddth0 = zeros(evalPoints,1);
P0 = [ss0; dss0; ddss0; th0; dth0; ddth0];

% objective function
fun = @(P) sum(1./P(numel(P)/6+1:2*numel(P)/6,1));

% linear inequality constraints
% monotonic path parameter
Aineq1 = diag(ones(evalPoints,1),0);
Aineq2 = diag(-1.*ones(evalPoints,1),1);
Aineq = Aineq1(1:evalPoints-1,1:evalPoints) + ...
    Aineq2(1:evalPoints-1,1:evalPoints);
Aineq = [Aineq zeros(size(Aineq,1),evalPoints*5)];
bineq = zeros(evalPoints-1,1);

% linear equality constraints
Aeq = zeros(6,length(P0));
Aeq(1,1) = 1; % position endpoints
Aeq(2,evalPoints) = 1;
Aeq(3,evalPoints+1) = fnval(dxp,sBounds(1)); % velocity endpoints
Aeq(4,2*evalPoints) = fnval(dxp,sBounds(2));
Aeq(5,evalPoints+1) = fnval(dyp,sBounds(1));
Aeq(6,2*evalPoints) = fnval(dyp,sBounds(2));
Aeq(7,3*evalPoints+1) = 1; % angle endpoints
Aeq(8,4*evalPoints) = 1;
Aeq(9,4*evalPoints+1) = 1; % angular velocity endpoints
Aeq(10,5*evalPoints) = 1;
beq = [sBounds(1); sBounds(2); ...
        s0(2); dx_des(1); ...
        s0(4); dx_des(2); ...
        s0(5); x_des(3); ...
        s0(6); dx_des(3)];

% lower and upper bounds
lb = [ones(evalPoints,1).*sBounds(1); ...
      zeros(evalPoints,1); ones(evalPoints,1).*-Inf; ...
      ones(evalPoints,1).*-pi; ones(evalPoints,1).*-Inf; ...
      ones(evalPoints,1).*-Inf];
ub = [ones(evalPoints,1).*sBounds(2); ...
      ones(evalPoints,1).*Inf; ones(evalPoints,1).*Inf; ...
      ones(evalPoints,1).*pi; ones(evalPoints,1).*Inf; ...
      ones(evalPoints,1).*Inf];

% nonlinear constraint function
nonlcon = @(P) nonlinconTOPP(P, s0, param, fCone, vec, tol, dxp, ddxp, dyp, ddyp);

% optimization
problem.x0 = P0;
problem.objective = fun;
problem.Aineq = Aineq;
problem.bineq = bineq;
problem.Aeq = Aeq;
problem.beq = beq;
problem.lb = lb;
problem.ub = ub;
problem.nonlcon = nonlcon;
problem.solver = 'fmincon';
problem.options = optimoptions('fmincon','Display','iter');

% invoke solver
psolve = fmincon(problem);

%% Results
P = reshape(psolve,numel(psolve)/6,6);