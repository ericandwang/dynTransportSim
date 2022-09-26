function [t0, t1] = intermediatePlan(th0, dth0, fCone, vec, type)
% Generating a trajectory that takes an statically infeasible IC to a
% statically feasible intermediate point (type 1) OR a statically feasible
% intermediate point to a statically infeasible FC (type -1)

if (sign(type) > 0)
    if (dth0 < 0 && th0 <= 1/2*dth0^2) || ...
            (dth0 >= 0 && th0 < -1/2*dth0^2)
        ddth0 = 1;
        ddth1 = -ddth0;
        t0 = roots([1/2*(ddth0^2+ddth0) dth0*ddth0 + dth0 th0 + 1/2*dth0^2]);
        t0 = max(t0);
        dth1 = dth0 + ddth0*t0;
        th1 = th0 + dth0*t0 + 1/2*ddth0*t0^2;
        t1 = abs(dth1/ddth1);
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
else
    t0 = 0;
    t1 = 0;
end

end

