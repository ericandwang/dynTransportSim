function DC = gradGenTOPP(r_GC, param, fCone, vec, s, dxp, ddxp, dyp, ddyp)
% NOTE: decomissioned

% parameters
g = param(9);
R_GC = norm(r_GC);
th_GC = angle(r_GC(1) + 1i*r_GC(2));

% symbolics
syms ds dds th dth ddth
ax_G = fnval(ddxp,s)*ds^2 + fnval(dxp,s)^2*dds;
ay_G = fnval(ddyp,s)*ds^2 + fnval(dyp,s)^2*dds;

p = [cos(th) sin(th) -R_GC*sin(th_GC); ...
              -sin(th) cos(th) R_GC*cos(th_GC); ...
              0 0 1]*[ax_G; ay_G; ddth] + ...
              [g*sin(th)-R_GC*dth^2*cos(th_GC); ...
              g*cos(th)-R_GC*dth^2*sin(th_GC); 0];

% calculating nonlinear constraint form
dim = size(fCone,2);
for i = 1:dim
    c(i,1) = vec(1,i)*(p(1)-fCone(1,i)) + vec(2,i)*(p(2)-fCone(2,i)) + ...
            vec(3,i)*(p(3)-fCone(3,i));
end

% calculating jacobian
states = [ds; dds; th; dth; ddth];
jac = jacobian(c,states)';

% generating function for jacobian
DC = matlabFunction(jac,'Vars',{states});

end