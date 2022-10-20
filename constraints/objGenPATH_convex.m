function [df] = objGenPATH_convex(P, r_GC, param, ss, knotVec, numCoefs, basisVectors, constraintSlopes)
% y distance from plane

% spline construction
syms coefs
coefs = sym('coefs',[numCoefs,1]);
numCoefs = length(coefs);
xcoefs = coefs(1:numCoefs/2);
ycoefs = coefs(numCoefs/2+1:end);
[~, ~, dxp_, ddxp_] = splineDiscreteReconstruct(knotVec, xcoefs, ss);
[~, ~, dyp_, ddyp_] = splineDiscreteReconstruct(knotVec, ycoefs, ss);

% constants
R_GC = norm(r_GC);
th_GC = angle(r_GC(1) + 1i*r_GC(2));
g = param(9);

% Objective function
f = 0;
for i = 1:length(ss)
    ds = P(i,1);
    dds = P(i,2);
    th = P(i,3);
    dth = P(i,4);
    ddth = P(i,5);
    ax_G = ddxp_(i)*ds^2 + dxp_(i)*dds;
    ay_G = ddyp_(i)*ds^2 + dyp_(i)*dds;
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

jac = jacobian(f,coefs)';

% generating function for objective gradient
df = matlabFunction(jac,'Vars',{coefs});

end

