function [xp, yp] = intermediatePlanStatic(xpD, ypD, knotVec, porder, direction)
% static path that connects to a dynamic path. Connects to end if direction
% = 1 and connect to start if direction = -1. Currently using straight line
% approach

xdisp = abs((fnval(xpD,1) - fnval(xpD,0))*1.5); % arbitrary distance
numPoints = 100;
dxpD = fnder(xpD,1);
dypD = fnder(ypD,1);

if (direction == 1)
    vAngRatio = fnval(dypD,1)/fnval(dxpD,1);
    xend = fnval(xpD,1) + xdisp*sign(fnval(dxpD,1));
    yend = fnval(ypD,1) + xdisp*sign(fnval(dxpD,1))*vAngRatio;
    x = linspace(fnval(xpD,1),xend,numPoints);
    y = linspace(fnval(ypD,1),yend,numPoints);
    s = linspace(0,1,length(x));
    xp = spap2(knotVec, porder, s, x);
    yp = spap2(knotVec, porder, s, y);
else
    vAngRatio = fnval(dypD,0)/fnval(dxpD,0);
    xend = fnval(xpD,0) - xdisp*sign(fnval(dxpD,0));
    yend = fnval(ypD,0) - xdisp*sign(fnval(dxpD,0))*vAngRatio;
    x = linspace(xend,fnval(xpD,0),numPoints);
    y = linspace(yend,fnval(ypD,0),numPoints);
    s = linspace(0,1,length(x));
    xp = spap2(knotVec, porder, s, x);
    yp = spap2(knotVec, porder, s, y);
end

% New integration method CCC TODO
if (direction == 1)
    ddxp = spap2(knotVec(3:end-2), porder-2, s, [zeros(1,99) -1]);
    dxp = fnint(ddxp,1);
    dxp.coefs = dxp.coefs - dxp.coefs(end);
    xp = fnint(dxp,0);
    xp.coefs = xp.coefs*xdisp/fnval(xp,1);
    yp = xp;
    yp.coefs = yp.coefs*vAngRatio;
    xp.coefs = xp.coefs + fnval(xpD,1);
    yp.coefs = yp.coefs + fnval(ypD,1);
else
    ddxp = spap2(knotVec(3:end-2), porder-2, s, [1 zeros(1,99)]);
    dxp = fnint(ddxp,1);
    dxp.coefs = dxp.coefs - dxp.coefs(1);
    xp = fnint(dxp,0);
    xp.coefs = xp.coefs*xdisp/fnval(xp,1);
    yp = xp;
    yp.coefs = yp.coefs*vAngRatio;
    xp.coefs = xp.coefs - xp.coefs(end) + fnval(xpD,1);
    yp.coefs = yp.coefs - yp.coefs(end) + fnval(ypD,1);
end

end

