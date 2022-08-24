function [c, ceq] = nonlinconPATH(P, r_GC, param, fCone, vec, ss, cc, endPoints, xxyy, porder, controlPoints, tol)

% control points spline construction
xx = [endPoints(1) xxyy(1:controlPoints-2) endPoints(2)];
yy = [endPoints(3) xxyy(controlPoints-2+1:end) endPoints(4)];
xp = spapi(optknt(cc,porder), cc, xx);
dxp = fnder(xp,1);
ddxp = fnder(xp,2);
yp = spapi(optknt(cc,porder), cc, yy);
dyp = fnder(yp,1);
ddyp = fnder(yp,2);

% constants
R_GC = norm(r_GC);
th_GC = angle(r_GC(1) + 1i*r_GC(2));
g = param(9);

% Gravito-inertial wrench
nPoints = length(ss);
for i = 1:nPoints
    ds = P(i,1);
    dds = P(i,2);
    th = P(i,3);
    dth = P(i,4);
    ddth = P(i,5);
    ax_G = fnval(ddxp,ss(i))*ds^2 + fnval(dxp,ss(i))*dds;
    ay_G = fnval(ddyp,ss(i))*ds^2 + fnval(dyp,ss(i))*dds;
    p(:,i) = [cos(th) sin(th) -R_GC*sin(th_GC); ...
              -sin(th) cos(th) R_GC*cos(th_GC); ...
              0 0 1]*[ax_G; ay_G; ddth] + ...
              [g*sin(th)-R_GC*dth^2*cos(th_GC); ...
              g*cos(th)-R_GC*dth^2*sin(th_GC); 0];
end

%% nonlinear inequality constraints
dim = size(fCone,2);
c = zeros(nPoints*dim,1);
for j = 1:nPoints
    for i = 1:dim
        c(i + dim*(j-1)) = dot(vec(:,i),p(:,j)-fCone(:,i)) + norm(vec(:,i))*tol;
    end
end


%% nonlinear equality constraints
ceq = [];

end
