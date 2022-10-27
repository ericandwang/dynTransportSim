function [th, dth, ddth, tTotal, xp, yp, x, vx, ax, y, vy, ay] = intermediatePlanDynamic(s0, x_des, fCone, param, knotVec, porder, type)
% Generating a trajectory that takes an statically infeasible IC to a
% statically feasible intermediate point (type 1) OR a statically feasible
% intermediate point to a statically infeasible FC (type -1)

% Parameters
g = param(9);
mu = param(10);

% Defining initial and desired values
x0 = s0(1);
dx0 = s0(2);
y0 = s0(3);
dy0 = s0(4);
th0 = s0(5);
dth0 = s0(6);
r_GC = s0(7:8,1);
xd = x_des(1);
yd = x_des(2);

% defining number of discrete points and linear acceleration limit and u (bang bang angular acceleration) CCC user input
u = 20;
numPoints = length(knotVec) - (porder-1)*2;
linAccelLim = 80;

% constants
R_GC = norm(r_GC);
th_GC = angle(r_GC(1) + 1i*r_GC(2));


% order of bang bang controller
accelOrder = 0;

if (dth0 < 0 && th0 <= 1/(2*u)*dth0^2) || ...
        (dth0 >= 0 && th0 < -1/(2*u)*dth0^2)
    accelOrder = 1;
    ddth0 = u;
    ddth1 = -ddth0;
    t0 = roots([ddth0, ...
                dth0*ddth0 + dth0, ...
                th0 + 1/(2*ddth0)*dth0^2]);
    t0 = max(t0);
    dth1 = dth0 + ddth0*t0;
    th1 = th0 + dth0*t0 + 1/2*ddth0*t0^2;
    t1 = t0 + abs(dth1/ddth1);
    tTotal = t1;
    th = @(t)  ((t >= 0) .* (t < t0)).*(th0 + dth0*t + 1/2*ddth0*t.^2) + ...
               ((t >= t0) .* (t < t1)).*(th1 + dth1*(t-t0) + 1/2*ddth1*(t-t0).^2) + ...
               (t >= t1).*(th1 + dth1*(t1-t0) + 1/2*ddth1*(t1-t0).^2);
    dth = @(t) ((t >= 0) .* (t < t0)).*(dth0 + ddth0*t) + ...
               ((t >= t0) .* (t < t1)).*(dth1 + ddth1*(t-t0)) + ...
               (t >= t1).*(dth1 + ddth1*(t1-t0));
    ddth = @(t) ((t >= 0) .* (t < t0)).*ddth0 + ...
                ((t >= t0) .* (t < t1)).*ddth1;
elseif (th0 == 0 && dth0 == 0)
    ddth0 = 0;
    ddth1 = 0;
    t0 = 0;
    t1 = 0;
else
    accelOrder = -1;
    ddth0 = -u;
    ddth1 = -ddth0;
    t0 = roots([1/2*(-ddth0^2+ddth0) -dth0*ddth0 + dth0 th0 - 1/2*dth0^2]);
    t0 = max(t0);
    dth1 = dth0 + ddth0*t0;
    th1 = th0 + dth0*t0 + 1/2*ddth0*t0^2;
    t1 = t0 + abs(dth1/ddth1);
    tTotal = t1;
    th = @(t)  ((t >= 0) .* (t < t0)).*(th0 + dth0*t + 1/2*ddth0*t.^2) + ...
               ((t >= t0) .* (t < t1)).*(th1 + dth1*(t-t0) + 1/2*ddth1*(t-t0).^2) + ...
               (t >= t1).*(th1 + dth1*(t1-t0) + 1/2*ddth1*(t1-t0).^2);
    dth = @(t) ((t >= 0) .* (t < t0)).*(dth0 + ddth0*t) + ...
               ((t >= t0) .* (t < t1)).*(dth1 + ddth1*(t-t0)) + ...
               (t >= t1).*(dth1 + ddth1*(t1-t0));
    ddth = @(t) ((t >= 0) .* (t < t0)).*ddth0 + ...
                ((t >= t0) .* (t < t1)).*ddth1;
end

% Discretizing time
t = linspace(0,tTotal,numPoints);
dt = t(2) - t(1);

% Defining boundary points for QP
boundaryPoints = zeros(2,4);
if (accelOrder == 1)
    boundaryPoints(:,1) = fCone(1:2,4)*u*accelOrder/fCone(3,4);
    boundaryPoints(:,2) = fCone(1:2,3)*u*accelOrder/fCone(3,3);
    boundaryPoints(:,3) = fCone(1:2,1)*u*-accelOrder/fCone(3,1);
    boundaryPoints(:,4) = fCone(1:2,2)*u*-accelOrder/fCone(3,2);
elseif (accelOrder == -1)
    boundaryPoints(:,1) = fCone(1:2,1)*u*accelOrder/fCone(3,1);
    boundaryPoints(:,2) = fCone(1:2,2)*u*accelOrder/fCone(3,2);
    boundaryPoints(:,3) = fCone(1:2,4)*u*-accelOrder/fCone(3,4);
    boundaryPoints(:,4) = fCone(1:2,3)*u*-accelOrder/fCone(3,3);
else
end
boundarySlopes = [(boundaryPoints(2,2)-boundaryPoints(2,1))/(boundaryPoints(1,2)-boundaryPoints(1,1)); ...
    (boundaryPoints(2,4)-boundaryPoints(2,3))/(boundaryPoints(1,4)-boundaryPoints(1,3));
    fCone(2,2)/fCone(1,2)];

% accelerations to gravito inertial point transformations
R_ = @(th) [cos(th) sin(th); -sin(th) cos(th)];
Rb_ = @(th,dth,ddth) ddth.*[-R_GC*sin(th_GC);R_GC*cos(th_GC)] + ...
    [g*sin(th)-R_GC*dth^2*cos(th_GC); g*cos(th)-R_GC*dth^2*sin(th_GC)];

% line constraints
Aline1 = [-boundarySlopes(3) -1; boundarySlopes(3) -1; boundarySlopes(1) -1];
bline1 = -[boundaryPoints(2,1)-(-boundarySlopes(3))*boundaryPoints(1,1); ...
          boundaryPoints(2,2)-boundarySlopes(3)*boundaryPoints(1,2); ...
          boundaryPoints(2,1)-boundarySlopes(1)*boundaryPoints(1,1)];
Aline2 = [-boundarySlopes(3) -1; boundarySlopes(3) -1; boundarySlopes(2) -1];
bline2 = -[boundaryPoints(2,3)-(-boundarySlopes(3))*boundaryPoints(1,3); ...
          boundaryPoints(2,4)-boundarySlopes(3)*boundaryPoints(1,4); ...
          boundaryPoints(2,1)-boundarySlopes(2)*boundaryPoints(1,1)];

% combining line constraints with transformations
A1 = @(th) Aline1*R_(th);
b1 = @(th,dth,ddth) bline1 - Aline1*Rb_(th,dth,ddth);
A2 = @(th) Aline2*R_(th);
b2 = @(th,dth,ddth) bline2 - Aline2*Rb_(th,dth,ddth);

% creating diagonal matrix and concatenating vector for constraints
Acell = cell(numPoints,1);
b = [];
for i = 1:numPoints
    if (t(i) < t0)
        Acell{i} = A1(th(t(i)));
        b = [b; b1(th(t(i)),dth(t(i)),ddth(t(i)))];
    else
        Acell{i} = A2(th(t(i)));
        b = [b; b2(th(t(i)),dth(t(i)),ddth(t(i)))];
    end
end
A = blkdiag(Acell{:});

prob = optimproblem('ObjectiveSense','min');
a = optimvar('a',numPoints*2,1,'LowerBound',ones(numPoints*2,1)*-linAccelLim, ...
    'UpperBound', ones(numPoints*2,1)*linAccelLim);
prob.Constraints.cons1 = A*a*type <= b;
objx = -xd + x0 + numPoints*dx0*dt;
objy = -yd + y0 + numPoints*dy0*dt;
for i = 1:numPoints
    objx = objx + (1/2+numPoints-i)*a(2*i-1)*dt^2;
    objy = objy + (1/2+numPoints-i)*a(2*i)*dt^2;
end
prob.Objective = objx^2 + objy^2;
opts = optimoptions('lsqlin');
%a0 = zeros(numPoints*2,1);
%a00 = struct('a',a0);
[asol, value] = solve(prob,'options',opts);

% Retrieving linear acceleration solution
ax = asol.a(1:2:end);
ay = asol.a(2:2:end);

% Plotting gravito inertial wrench constraints
figure(200)
hold off
patch([0; fCone(1,1); fCone(1,2)],[0;fCone(2,1);fCone(2,2)],[0;fCone(3,1);fCone(3,2)],[0.5,0,0.5],'FaceAlpha',0.1)
hold on
patch([0; fCone(1,2); fCone(1,3)],[0;fCone(2,2);fCone(2,3)],[0;fCone(3,2);fCone(3,3)],[0.5,0,0.5],'FaceAlpha',0.1)
patch([0; fCone(1,3); fCone(1,4)],[0;fCone(2,3);fCone(2,4)],[0;fCone(3,3);fCone(3,4)],[0.5,0,0.5],'FaceAlpha',0.1)
patch([0; fCone(1,4); fCone(1,1)],[0;fCone(2,4);fCone(2,1)],[0;fCone(3,4);fCone(3,1)],[0.5,0,0.5],'FaceAlpha',0.1)
p = zeros(3,numPoints-1);
for i = 1:numPoints-1
    th1 = th(t(i));
    dth1 = dth(t(i));
    ddth1 = ddth(t(i));
    ax_G = ax(i);
    ay_G = ay(i);
    p(:,i) = [cos(th1) sin(th1) -R_GC*sin(th_GC); ...
              -sin(th1) cos(th1) R_GC*cos(th_GC); ...
              0 0 1]*[ax_G; ay_G; ddth1] + ...
              [g*sin(th1)-R_GC*dth1^2*cos(th_GC); ...
              g*cos(th1)-R_GC*dth1^2*sin(th_GC); 0];
end
scatter3(p(1,1),p(2,1),p(3,1),'k');
plot3(p(1,:),p(2,:),p(3,:),'b')
xlabel('$\ddot{x}_{c} [m/s^2]$','interpreter','latex')
ylabel('$\ddot{y}_{c} [m/s^2]$','interpreter','latex')
zlabel('$\ddot{\theta}_{c} [rad/s^2]$','interpreter','latex')
title('Gravito-Inertial Wrench Constraints')

x = zeros(numPoints,1);
vx = zeros(numPoints,1);
x(1) = x0;
vx(1) = dx0;
y = zeros(numPoints,1);
vy = zeros(numPoints,1);
y(1) = y0;
vy(1) = dy0;
for i = 1:length(x)-1
    x(i+1) = x(i) + vx(i)*dt + 1/6*(2*ax(i)+ax(i+1))/2*dt^2;
    vx(i+1) = vx(i) + (ax(i)+ax(i+1))/2*dt;
    y(i+1) = y(i) + vy(i)*dt + 1/6*(2*ay(i)+ay(i+1))*dt^2;
    vy(i+1) = vy(i) + (ay(i)+ay(i+1))*dt;
end

% Plotting path debug
%figure(1001), hold off, plot(x,y);
%hold on, scatter(x(floor(length(x)/2)),y(floor(length(y)/2)))

% Obtaining path
s = linspace(0,1,length(x));
%xp = spap2(knotVec, porder, s, x);
%yp = spap2(knotVec, porder, s, y);
ddxp = spap2(knotVec(3:end-2), porder-2, s, ax*tTotal^2);
ddyp = spap2(knotVec(3:end-2), porder-2, s, ay*tTotal^2);
dxp = fnint(ddxp,vx(1)*tTotal);
dyp = fnint(ddyp,vy(1)*tTotal);
xp = fnint(dxp,x(1));
yp = fnint(dyp,y(1));

if (type == -1) % flipping representation
    xp.coefs = flip(xp.coefs);
    yp.coefs = flip(yp.coefs);
    x = flip(x); vx = flip(vx); ax = flip(ax);
    y = flip(y); vy = flip(vy); ay = flip(ay);
    th = @(t) th(tTotal - t);
    dth = @(t) dth(tTotal - t);
    ddth = @(t) ddth(tTotal - t);
end

% Appending path to drive vx and vy to 0 (statically stable section)
% CCC causes unbounded path to exit workspace due to limited
% acceleration UNUSED DELETE
%     xend = x(end);
%     yend = y(end);
%     vxend = vx(end);
%     vyend = vy(end);
%     if (vxend > 0 && vyend > 0)
%         A = [-1 -mu; -1 vxend/vyend];
%         b = [mu*g; 0];
%         a1 = A\b;
%     else
%         a1 = [0;0]; % CCC fix for other cases
%     end
%     ax1 = a1(1); ay1 = a1(2);
%     tTotal1 = abs(vxend/ax1);
%     x1end = xend + vxend*tTotal1 + 1/2*ax1*tTotal1^2;
%     y1end = yend + vyend*tTotal1 + 1/2*ay1*tTotal1^2;



end

