function [th, dth, ddth, tTotal] = intermediatePlan(th0, dth0, fCone, vec, r_GC, g, type)
% Generating a trajectory that takes an statically infeasible IC to a
% statically feasible intermediate point (type 1) OR a statically feasible
% intermediate point to a statically infeasible FC (type -1)

% defining magnitude of angular acceleration CCC turn into parameter
angAccel = 1;
numPoints = 100;

% constants
R_GC = norm(r_GC);
th_GC = angle(r_GC(1) + 1i*r_GC(2));


% order of bang bang controller
accelOrder = 0;

if (sign(type) > 0) % going from IC to intermediate point
    if (dth0 < 0 && th0 <= 1/2*dth0^2) || ...
            (dth0 >= 0 && th0 < -1/2*dth0^2)
        accelOrder = 1;
        ddth0 = angAccel;
        ddth1 = -ddth0;
        t0 = roots([1/2*(ddth0^2+ddth0) dth0*ddth0 + dth0 th0 + 1/2*dth0^2]);
        t0 = max(t0);
        dth1 = dth0 + ddth0*t0;
        th1 = th0 + dth0*t0 + 1/2*ddth0*t0^2;
        t1 = t0 + abs(dth1/ddth1);
        tTotal = t1;
        th = @(t)  ((t >= 0) .* (t < t0)).*(th0 + dth0*t + 1/2*ddth0*t.^2) + ...
                   ((t >= t0) .* (t < t1)).*(th1 + dth1*(t-t0) + 1/2*ddth1*(t-t0).^2) + ...
                   (t >= t1).*(th1 + dth1*(t1-t0) + 1/2*ddth1*(t1-t0).^2);
        dth = @(t) ((t >= 0) .* (t < t0)).*(dth0 + ddth0*t) + ...
                   ((t >= t0) .* (t < t1)).*(dth1 + ddth1*(t-t0)) + ...
                   (t >= t1).*(dth1 + ddth1*(t1-t0));
        ddth = @(t) ((t >= 0) .* (t < t0)).*ddth0 + ...
                    ((t >= t0) .* (t < t1)).*ddth1;
    elseif (th0 == 0 && dth0 == 0)
        ddth0 = 0;
        ddth1 = 0;
        t0 = 0;
        t1 = 0;
    else
        accelOrder = -1;
        ddth0 = -angAccel;
        ddth1 = -ddth0;
        t0 = roots([1/2*(-ddth0^2+ddth0) -dth0*ddth0 + dth0 th0 - 1/2*dth0^2]);
        t0 = max(t0);
        dth1 = dth0 + ddth0*t0;
        th1 = th0 + dth0*t0 + 1/2*ddth0*t0^2;
        t1 = t0 + abs(dth1/ddth1);
        tTotal = t1;
        th = @(t)  ((t >= 0) .* (t < t0)).*(th0 + dth0*t + 1/2*ddth0*t.^2) + ...
                   ((t >= t0) .* (t < t1)).*(th1 + dth1*(t-t0) + 1/2*ddth1*(t-t0).^2) + ...
                   (t >= t1).*(th1 + dth1*(t1-t0) + 1/2*ddth1*(t1-t0).^2);
        dth = @(t) ((t >= 0) .* (t < t0)).*(dth0 + ddth0*t) + ...
                   ((t >= t0) .* (t < t1)).*(dth1 + ddth1*(t-t0)) + ...
                   (t >= t1).*(dth1 + ddth1*(t1-t0));
        ddth = @(t) ((t >= 0) .* (t < t0)).*ddth0 + ...
                    ((t >= t0) .* (t < t1)).*ddth1;
    end

    % Discretizing time
    t = linspace(0,tTotal,numPoints);

    % Defining boundary points for QP
    boundaryPoints = zeros(2,4);
    if (accelOrder == 1)
        boundaryPoints(:,1) = fCone(1:2,4)*angAccel*accelOrder/fCone(3,4);
        boundaryPoints(:,2) = fCone(1:2,3)*angAccel*accelOrder/fCone(3,3);
        boundaryPoints(:,3) = fCone(1:2,1)*angAccel*accelOrder/fCone(3,1);
        boundaryPoints(:,4) = fCone(1:2,2)*angAccel*accelOrder/fCone(3,2);
    elseif (accelOrder == -1)
        boundaryPoints(:,1) = fCone(1:2,1)*angAccel*accelOrder/fCone(3,1);
        boundaryPoints(:,2) = fCone(1:2,2)*angAccel*accelOrder/fCone(3,2);
        boundaryPoints(:,3) = fCone(1:2,4)*angAccel*accelOrder/fCone(3,4);
        boundaryPoints(:,4) = fCone(1:2,3)*angAccel*accelOrder/fCone(3,3);
    else
    end
    boundarySlopes = [(boundaryPoints(2,2)-boundaryPoints(2,1))/(boundaryPoints(1,2)-boundaryPoints(1,1)); ...
        (boundaryPoints(2,4)-boundaryPoints(2,3))/(boundaryPoints(1,4)-boundaryPoints(1,3));
        fCone(2,2)/fCone(1,2)];

    % accelerations to gravito inertial point transformations
    R_ = @(th) [cos(th) sin(th); -sin(th) cos(th)];
    Rb_ = @(th,dth,ddth) ddth.*[-R_GC*sin(th_GC);R_GC*cos(th_GC)] + ...
        [g*sin(th)-R_GC*dth^2*cos(th_GC); g*cos(th)-R_GC*dth^2*sin(th_GC)];

    % line constraints
    Aline1 = [-boundarySlopes(3) -1; boundarySlopes(3) -1; boundarySlopes(1) -1];
    bline1 = [boundaryPoints(2,1)-(-boundarySlopes(3))*boundaryPoints(1,1); ...
              boundaryPoints(2,2)-boundarySlopes(3)*boundaryPoints(1,2); ...
              boundaryPoints(2,1)-boundarySlopes(1)*boundaryPoints(1,1)];
    Aline2 = [-boundarySlopes(3) -1; boundarySlopes(3) -1; boundarySlopes(2) -1];
    bline2 = [boundaryPoints(2,3)-(-boundarySlopes(3))*boundaryPoints(1,3); ...
              boundaryPoints(2,4)-boundarySlopes(3)*boundaryPoints(1,4); ...
              boundaryPoints(2,1)-boundarySlopes(2)*boundaryPoints(1,1)];

    % combining line constraints with transformations
    A1 = @(th) Aline1*R_(th);
    b1 = @(th,dth,ddth) bline1 - Aline1*Rb_(th,dth,ddth);
    A2 = @(th) Aline2*R_(th);
    b2 = @(th,dth,ddth) bline2 - Aline2*Rb_(th,dth,ddth);

    % creating diagonal matrix and concatenating vector for constraints
    Acell = cell(numPoints,1);
    b = [];
    for i = 1:numPoints
        if (t(i) < t0)
            Acell{i} = A1(th(t(i)));
            b = [b; b1(th(t(i)),dth(t(i)),ddth(t(i)))];
        else
            Acell{i} = A2(th(t(i)));
            b = [b; b2(th(t(i)),dth(t(i)),ddth(t(i)))];
        end
    end
    A = blkdiag(Acell{:});


    % CCC change to QP determining discrete ax_G and ay_G trajectory
%     p(:,ii) = [cos(th) sin(th) -R_GC*sin(th_GC); ...
%               -sin(th) cos(th) R_GC*cos(th_GC); ...
%               0 0 1]*[ax_G; ay_G; ddth] + ...
%               [g*sin(th)-R_GC*dth^2*cos(th_GC); ...
%               g*cos(th)-R_GC*dth^2*sin(th_GC); 0];
%     nPoints = size(p,2);
%     dim = size(fCone,2);
%     c = zeros(nPoints*dim,1);
%     for j = 1:nPoints
%         for i = 1:dim
%             c(i + dim*(j-1)) = dot(vec(:,i),p(:,j)-fCone(:,i)) + norm(vec(:,i))*tol;
%         end
%     end

else % CCC TODO going from intermediate point to FC
    t0 = 0;
    t1 = 0;
end


%% Helper functions

    function hi = hello(x)
        hi = x^2;
    end

end

