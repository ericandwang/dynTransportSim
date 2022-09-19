function [c, ceq, dc, dceq] = nonlinconTOPP(P, s0, ss, param, fCone, vec, tol, xp, dxp, ddxp, yp, dyp, ddyp, dcFun, dceqFun, accelLim)

% initial states
r_GC = [s0(7); s0(8)]; % object frame

% parameters
g = param(9);

% states
Pmat = P;
P = reshape(P,numel(P)/5,5);
dss = P(:,1);
ddss = P(:,2);
ths = P(:,3);
dths = P(:,4);
ddths = P(:,5);

%% nonlinear inequality constraints

% converting to gravito-inertial wrench pgi
nPoints = size(P,1);
R_GC = norm(r_GC);
th_GC = angle(r_GC(1) + 1i*r_GC(2));
xp_ = fnval(xp,ss);
dxp_ = fnval(dxp,ss);
ddxp_ = fnval(ddxp,ss);
yp_ = fnval(yp,ss);
dyp_ = fnval(dyp,ss);
ddyp_ = fnval(ddyp,ss);
for ii = 1:nPoints
    ds = dss(ii);
    dds = ddss(ii);
    th = ths(ii);
    dth = dths(ii);
    ddth = ddths(ii);
    ax_G = ddxp_(ii)*ds^2 + dxp_(ii)*dds;
    ay_G = ddyp_(ii)*ds^2 + dyp_(ii)*dds;
    p(:,ii) = [cos(th) sin(th) -R_GC*sin(th_GC); ...
              -sin(th) cos(th) R_GC*cos(th_GC); ...
              0 0 1]*[ax_G; ay_G; ddth] + ...
              [g*sin(th)-R_GC*dth^2*cos(th_GC); ...
              g*cos(th)-R_GC*dth^2*sin(th_GC); 0];
end

% appending collocation points
h = ss(2)-ss(1);
for ii = 1:nPoints-1
    dds = ddss(ii);
    ds = sqrt(h*dds + dss(ii)^2);
    dx = -3/(2*h)*(xp_(ii)-xp_(ii+1)) - 1/4*(dxp_(ii)+dxp_(ii+1));
    dy = -3/(2*h)*(yp_(ii)-yp_(ii+1)) - 1/4*(dyp_(ii)+dyp_(ii+1));
    ddx = (ddxp_(ii)+ddxp_(ii+1))/2;
    ddy = (ddyp_(ii)+ddyp_(ii+1))/2;
    dt = h/(dss(ii)+ds);
    ddth = ddths(ii);
    dth = dths(ii) + ddth*dt;
    th = ths(ii) + dth*dt + 1/2*ddth*dt^2;
    ax_G = ddx*ds^2 + dx*dds;
    ay_G = ddy*ds^2 + dy*dds;
    pappend(:,ii) = [cos(th) sin(th) -R_GC*sin(th_GC); ...
              -sin(th) cos(th) R_GC*cos(th_GC); ...
              0 0 1]*[ax_G; ay_G; ddth] + ...
              [g*sin(th)-R_GC*dth^2*cos(th_GC); ...
              g*cos(th)-R_GC*dth^2*sin(th_GC); 0];
end

p = [p, pappend];


% constructing linear plane constraints (knot points and collocation points)
nPoints = size(p,2);
dim = size(fCone,2);
c = zeros(nPoints*dim,1);
for j = 1:nPoints
    for i = 1:dim
        c(i + dim*(j-1)) = dot(vec(:,i),p(:,j)-fCone(:,i)) + norm(vec(:,i))*tol;
    end
end

% appending object frame y acceleration
cappend = p(2,:)' - accelLim;
c = [c; cappend];

%% nonlinear equality constraints

% path acceleration constraints
ceq1 = dss(2:end).^2 - dss(1:end-1).^2 - ...
    2.*(ss(2:end)-ss(1:end-1)).*ddss(1:end-1);

% orientation kinematic constraints
dt = 2.*(ss(2:end)-ss(1:end-1))./(dss(2:end)+dss(1:end-1));
ceq2 = dths(2:end)-dths(1:end-1) - ddths(1:end-1).*dt;
ceq3 = (ths(1:end-1)-ths(2:end)) + dths(1:end-1).*dt + ...
    1/2.*ddths(1:end-1).*dt.^2;

ceq = [ceq1; ceq2; ceq3];

%% gradients
if nargout > 2
dc = dcFun(Pmat);

% equality gradients
dceq = dceqFun(Pmat);
end

end

