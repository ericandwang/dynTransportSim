function [c, ceq] = nonlinconTOPP(P, s0, ss0, param, fCone, vec, tol, dxp, ddxp, dyp, ddyp)

% initial states
r_GC = [s0(7); s0(8)]; % object frame

% parameters
g = param(9);

% states
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
    th = ths(ii);
    dth = dths(ii);
    ddth = ddths(ii);
    s = ss(ii);
    ds = dss(ii);
    dds = ddss(ii);
    ax_G = fnval(ddxp,s)*ds^2 + fnval(dxp,s)^2*dds;
    ay_G = fnval(ddyp,s)*ds^2 + fnval(dyp,s)^2*dds;
    p(:,ii) = [cos(th) sin(th) -R_GC*sin(th_GC); ...
              -sin(th) cos(th) R_GC*cos(th_GC); ...
              0 0 1]*[ax_G; ay_G; ddth] + ...
              [g*sin(th)-R_GC*dth^2*cos(th_GC); ...
              g*cos(th)-R_GC*dth^2*sin(th_GC); 0];
end


% constructing linear plane constraints
dim = size(fCone,2);
c = zeros(nPoints*dim,1);
for j = 1:nPoints
    for i = 1:dim
        c(i + dim*(j-1)) = vec(1,i)*(p(1,j)-fCone(1,i)) + vec(2,i)*(p(2,j)-fCone(2,i)) + ...
            vec(3,i)*(p(3,j)-fCone(3,i)) + norm(vec(:,i))*tol;
    end
end

%% nonlinear equality constraints

% path acceleration constraints
ceq1 = dss(2:end).^2 - dss(1:end-1).^2 - ...
    2.*(ss(2:end)-ss(1:end-1)).*ddss(1:end-1);

% orientation kinematic constraints
%ceq2 = dths(2:end).^2 - dths(1:end-1).^2 - ...
%    2.*(ths(2:end)-ths(1:end-1)).*ddths(1:end-1);
dt = 2.*(ss(2:end)-ss(1:end-1))./(dss(2:end)+dss(1:end-1));
ceq2 = dths(2:end)-dths(1:end-1) - ddths(1:end-1).*dt;
ceq3 = (ths(1:end-1)-ths(2:end)) + dths(1:end-1).*dt + ...
    1/2.*ddths(1:end-1).*dt.^2;

ceq = [ceq1; ceq2; ceq3];

end

