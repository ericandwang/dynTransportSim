function [xp, yp] = intermediatePlanBridge(x0, y0, xf, yf, knotVec, porder)
% static path that connects to static paths on both ends.

numPoints = 100;
s = linspace(0,1,numPoints);
ddxp = spap2(knotVec(3:end-2), porder-2, s, [1 zeros(1,98) -1]);
dxp = fnint(ddxp,0);
xp = fnint(dxp,0);
xp.coefs = xp.coefs*(xf-x0)/fnval(xp,1);
yp = xp;
yp.coefs = yp.coefs*(yf-y0)/fnval(yp,1);
xp.coefs = xp.coefs + x0;
yp.coefs = yp.coefs + y0;

end

