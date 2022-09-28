function [th, dth, ddth, tTotal] = intermediatePlan(th0, dth0, fCone, vec, type)
% Generating a trajectory that takes an statically infeasible IC to a
% statically feasible intermediate point (type 1) OR a statically feasible
% intermediate point to a statically infeasible FC (type -1)

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
        t1 = abs(dth1/ddth1);
        t1 = t1;
    end
else % going from intermediate point to FC
    t0 = 0;
    t1 = 0;
end

end

