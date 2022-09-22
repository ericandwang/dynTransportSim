function dc = dcGenPATH(P, r_GC, param, fCone, vec, ss, knotVec, numCoefs, accelLim)

% spline construction
syms coefs
coefs = sym('coefs',[numCoefs,1]);
numCoefs = length(coefs);
xcoefs = coefs(1:numCoefs/2);
ycoefs = coefs(numCoefs/2+1:end);
[~, xp_, dxp_, ddxp_] = splineDiscreteReconstruct(knotVec, xcoefs);
[~, yp_, dyp_, ddyp_] = splineDiscreteReconstruct(knotVec, ycoefs);

% constants
R_GC = norm(r_GC);
th_GC = angle(r_GC(1) + 1i*r_GC(2));
g = param(9);

% Gravito-inertial wrench
nPoints = length(ss);
dss = P(:,1);
ddss = P(:,2);
ths = P(:,3);
dths = P(:,4);
ddths = P(:,5);
for i = 1:nPoints
    ds = dss(i);
    dds = ddss(i);
    th = ths(i);
    dth = dths(i);
    ddth = ddths(i);
    ax_G = ddxp_(i)*ds^2 + dxp_(i)*dds;
    ay_G = ddyp_(i)*ds^2 + dyp_(i)*dds;
    p(:,i) = [cos(th) sin(th) -R_GC*sin(th_GC); ...
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

%p = [p, pappend];

% constructing linear plane constraints
nPoints = size(p,2);
dim = size(fCone,2);
for j = 1:nPoints
    for i = 1:dim
        c(i + dim*(j-1),1) = dot(vec(:,i),p(:,j)-fCone(:,i));
    end
end

% appending object frame y acceleration
cappend = p(2,:)' - accelLim;
c = [c; cappend];


jac = jacobian(c,coefs)';

% generating function for inequality gradient
dc = matlabFunction(jac,'Vars',{coefs});

end

