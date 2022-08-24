function [c, ceq, dc, dceq] = nonlinconTOPP_lin(P, s0, ss0, param, fCone, vec, tol, dxp, ddxp, dyp, ddyp, dcFun, dceqFun, th0)

% initial states
r_GC = [s0(7); s0(8)]; % object frame

% parameters
g = param(9);

% states
Pmat = P;
P = reshape(P,numel(P)/5,5);
ss = ss0;
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
p = zeros(3,nPoints);
for ii = 1:nPoints
    s = ss(ii);
    ds = dss(ii);
    dds = ddss(ii);
    th = ths(ii);
    dth = dths(ii);
    ddth = ddths(ii);
    ax_G = fnval(ddxp,s)*ds^2 + fnval(dxp,s)*dds;
    ay_G = fnval(ddyp,s)*ds^2 + fnval(dyp,s)*dds;
    p(:,ii) = [cos(th0(ii))-sin(th0(ii))*(th-th0(ii)) sin(th0(ii))+cos(th0(ii))*(th-th0(ii)) -R_GC*sin(th_GC); ...
              -(sin(th0(ii))+cos(th0(ii))*(th-th0(ii))) cos(th0(ii))-sin(th0(ii))*(th-th0(ii)) R_GC*cos(th_GC); ...
              0 0 1]*[ax_G; ay_G; ddth] + ...
              [g*(sin(th0(ii))+cos(th0(ii))*(th-th0(ii)))-R_GC*dth^2*cos(th_GC); ...
              g*(cos(th0(ii))-sin(th0(ii))*(th-th0(ii)))-R_GC*dth^2*sin(th_GC); 0];
end


% constructing linear plane constraints
dim = size(fCone,2);
c = zeros(nPoints*dim,1);
for j = 1:nPoints
    for i = 1:dim
        c(i + dim*(j-1)) = dot(vec(:,i),p(:,j)-fCone(:,i)) + norm(vec(:,i))*tol;
    end
end

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

