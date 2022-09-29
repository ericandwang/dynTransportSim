function [th, dth, ddth, tTotal] = intermediatePlan(th0, dth0, fCone, vec, r_GC, type)
% Generating a trajectory that takes an statically infeasible IC to a
% statically feasible intermediate point (type 1) OR a statically feasible
% intermediate point to a statically infeasible FC (type -1)

% constants
R_GC = norm(r_GC);
th_GC = angle(r_GC(1) + 1i*r_GC(2));

if (sign(type) > 0) % going from IC to intermediate point
    if (dth0 < 0 && th0 <= 1/2*dth0^2) || ...
            (dth0 >= 0 && th0 < -1/2*dth0^2)
        ddth0 = 1;
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
        ddth0 = -1;
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

    % CCC change to QP determining discrete ax_G and ay_G trajectory
    p(:,ii) = [cos(th) sin(th) -R_GC*sin(th_GC); ...
              -sin(th) cos(th) R_GC*cos(th_GC); ...
              0 0 1]*[ax_G; ay_G; ddth] + ...
              [g*sin(th)-R_GC*dth^2*cos(th_GC); ...
              g*cos(th)-R_GC*dth^2*sin(th_GC); 0];
    nPoints = size(p,2);
    dim = size(fCone,2);
    c = zeros(nPoints*dim,1);
    for j = 1:nPoints
        for i = 1:dim
            c(i + dim*(j-1)) = dot(vec(:,i),p(:,j)-fCone(:,i)) + norm(vec(:,i))*tol;
        end
    end

else % CCC TODO going from intermediate point to FC
    t0 = 0;
    t1 = 0;
end

end

