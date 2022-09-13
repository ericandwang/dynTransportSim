function [f] = objPATH_convex(P, r_GC, param, fCone, vec, ss, knotVec, coefs, basisVectors, constraintSlopes)
% y distance from plane

% spline construction
numCoefs = length(coefs);
xcoefs = coefs(1:numCoefs/2);
ycoefs = coefs(numCoefs/2+1:end);
xp = spmak(knotVec,xcoefs');
dxp = fnder(xp,1);
ddxp = fnder(xp,2);
yp = spmak(knotVec,ycoefs');
dyp = fnder(yp,1);
ddyp = fnder(yp,2);

% constants
R_GC = norm(r_GC);
th_GC = angle(r_GC(1) + 1i*r_GC(2));
g = param(9);

% Objective function
f = 0;
dim = size(fCone,2);
for i = 1:length(ss)
    ds = P(i,1);
    dds = P(i,2);
    th = P(i,3);
    dth = P(i,4);
    ddth = P(i,5);
    ax_G = fnval(ddxp,ss(i))*ds^2 + fnval(dxp,ss(i))*dds;
    ay_G = fnval(ddyp,ss(i))*ds^2 + fnval(dyp,ss(i))*dds;
    p = [cos(th) sin(th) -R_GC*sin(th_GC); ...
              -sin(th) cos(th) R_GC*cos(th_GC); ...
              0 0 1]*[ax_G; ay_G; ddth] + ...
              [g*sin(th)-R_GC*dth^2*cos(th_GC); ...
              g*cos(th)-R_GC*dth^2*sin(th_GC); 0];
    p_ = basisVectors \ p;
    distx = p(2) - constraintSlopes(1)*abs(p(1));
    distth = p_(2) - constraintSlopes(2)*abs(p_(3));
    f = f - distx - distth;
end

end

