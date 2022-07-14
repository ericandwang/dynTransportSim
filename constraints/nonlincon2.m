function [c, ceq] = nonlincon2(nHorizon, x, s, param, fCone, vec, tol)

% initial states
r_GC = [s(7); s(8)]; % object frame

% parameters
g = param(9);

% converting to object accelerations p
R_GC = norm(r_GC);
th_GC = angle(r_GC(1) + 1i*r_GC(2));
p = zeros(3,nHorizon);
for i = 1:nHorizon
    th = x(5 + 9*(i-1));
    dth = x(6 + 9*(i-1));
    ax_G = x(7 + 9*(i-1));
    ay_G = x(8 + 9*(i-1));
    ath_G = x(9 + 9*(i-1));
    p(:,i) = [cos(th) sin(th) -R_GC*sin(th_GC); ...
              -sin(th) cos(th) R_GC*cos(th_GC); ...
              0 0 1]*[ax_G; ay_G; ath_G] + ...
              [g*sin(th)-R_GC*dth^2*cos(th_GC); ...
              g*cos(th)-R_GC*dth^2*sin(th_GC); 0];
end


% constructing linear plane constraints
% tolerance for distance from plane
%tol = 0.5; %0.1
dim = size(fCone,2);
c = zeros(nHorizon*dim,1);
for j = 1:nHorizon
    for i = 1:dim
        c(i + dim*(j-1)) = vec(1,i)*(p(1,j)-fCone(1,i)) + vec(2,i)*(p(2,j)-fCone(2,i)) + ...
            vec(3,i)*(p(3,j)-fCone(3,i)) + norm(vec(:,i))*tol;
    end
end

% no nonlinear equality constraints
ceq = [];

end

