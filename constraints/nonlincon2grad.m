function [c, ceq, DC, DCeq] = nonlincon2grad(nHorizon, x, s, param, fCone, vec, tol, jacob)

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

% deriving nonlinear inequality and equality constraints gradients
if nargout > 2
    DC = zeros(nHorizon*9,nHorizon*dim);
    for i = 1:nHorizon
        stateVec = x(1+9*(i-1):9+9*(i-1));
        DC(1+9*(i-1):9+9*(i-1),1+4*(i-1):4+4*(i-1)) = ...
            jacob(stateVec(1),stateVec(2),stateVec(3),stateVec(4), ...
            stateVec(5),stateVec(6),stateVec(7),stateVec(8),stateVec(9));
    end
    % no nonlinear equality constraints
    DCeq = [];
end

end

