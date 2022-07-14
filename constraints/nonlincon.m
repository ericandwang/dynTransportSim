function [c, ceq] = nonlincon(x, s, param, fCone, vec)

% states
th = s(5);
dth = s(6);
r_GC = [s(7); s(8)]; % object frame

% world accelerations
ax_G = x(1);
ay_G = x(2);
ath_G = x(3);

% parameters
g = param(9);

% converting to object acceleration p
R_GC = norm(r_GC);
th_GC = angle(r_GC(1) + 1i*r_GC(2));
p = [cos(th) sin(th) -R_GC*sin(th_GC); ...
 -sin(th) cos(th) R_GC*cos(th_GC); ...
 0 0 1]*[ax_G; ay_G; ath_G] + ...
    [g*sin(th)-R_GC*dth^2*cos(th_GC); ...
     g*cos(th)-R_GC*dth^2*sin(th_GC); 0];

% constructing linear plane constraints
% tolerance for distance from plane
tol = 1;
c = zeros(size(fCone,2),1);
for i = 1:size(fCone,2)
    c(i) = vec(1,i)*(p(1)-fCone(1,i)) + vec(2,i)*(p(2)-fCone(2,i)) + ...
        vec(3,i)*(p(3)-fCone(3,i)) + norm(vec(:,i))*tol;
end
ceq = [];

end

