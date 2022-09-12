function dc = dcGenPATH(P, r_GC, param, fCone, vec, ss, knotVec, numCoefs, accelLim)

% spline construction
syms coefs
coefs = sym('coefs',[numCoefs,1]);
numCoefs = length(coefs);
xcoefs = coefs(1:numCoefs/2);
ycoefs = coefs(numCoefs/2+1:end);
[~, ~, dxp_, ddxp_] = splineDiscreteReconstruct(knotVec, xcoefs);
[~, ~, dyp_, ddyp_] = splineDiscreteReconstruct(knotVec, ycoefs);

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
    ax_G = ddxp_(i)*ds^2 + dxp_(i)*dds;
    ay_G = ddyp_(i)*ds^2 + dyp_(i)*dds;
    p(:,i) = [cos(th) sin(th) -R_GC*sin(th_GC); ...
              -sin(th) cos(th) R_GC*cos(th_GC); ...
              0 0 1]*[ax_G; ay_G; ddth] + ...
              [g*sin(th)-R_GC*dth^2*cos(th_GC); ...
              g*cos(th)-R_GC*dth^2*sin(th_GC); 0];
end

% constructing linear plane constraints
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

