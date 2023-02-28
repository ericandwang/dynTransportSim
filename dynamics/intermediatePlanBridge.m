%function [xp, yp] = intermediatePlanBridge(x0, y0, xf, yf, knotVec, porder)
function [xp_2, yp_2] = intermediatePlanBridge(xp_1, yp_1, xp_3, yp_3, knotVec, porder)
% static path that connects to static paths on both ends.

% numPoints = 100;
% s = linspace(0,1,numPoints);
% ddxp = spap2(knotVec(3:end-2), porder-2, s, [1 zeros(1,98) -1]);
% dxp = fnint(ddxp,0);
% xp = fnint(dxp,0);
% xp.coefs = xp.coefs*(xf-x0)/fnval(xp,1);
% yp = xp;
% yp.coefs = yp.coefs*(yf-y0)/fnval(yp,1);
% xp.coefs = xp.coefs + x0;
% yp.coefs = yp.coefs + y0;


dxp_1 = fnder(xp_1,1);
dxp_3 = fnder(xp_3,1);
dyp_1 = fnder(yp_1,1);
dyp_3 = fnder(yp_3,1);
ddxp_1 = fnder(xp_1,2);
ddxp_3 = fnder(xp_3,2);
ddyp_1 = fnder(yp_1,2);
ddyp_3 = fnder(yp_3,2);

% assuming p_2 is same spline structure as p_1 and xp and yp both have
% identical spline knot structures
s1 = xp_1.knots(end);
s3 = xp_3.knots(1);
xp_2 = xp_1;
yp_2 = yp_1;
xp_2.knots = xp_2.knots + - xp_2.knots(1);
yp_2.knots = yp_2.knots + - yp_2.knots(1);
[~, ~, ~, posEndpoints, velEndpoints, accelEndpoints] = splineBasisCoefs(xp_1.knots, xp_1.order);
A1 = [posEndpoints(1,1:3); velEndpoints(1,1:3); accelEndpoints(1,1:3)];
A3 = [posEndpoints(2,end-2:end); velEndpoints(2,end-2:end); accelEndpoints(2,end-2:end)];
bx1 = [fnval(xp_1,s1); fnval(dxp_1,s1); fnval(ddxp_1,s1)];
bx3 = [fnval(xp_3,s3); fnval(dxp_3,s3); fnval(ddxp_3,s3)];
by1 = [fnval(yp_1,s1); fnval(dyp_1,s1); fnval(ddyp_1,s1)];
by3 = [fnval(yp_3,s3); fnval(dyp_3,s3); fnval(ddyp_3,s3)];
xp_2.coefs(1:3) = (A1 \ bx1)';
xp_2.coefs(end-2:end) = (A3 \ bx3)';
yp_2.coefs(1:3) = (A1 \ by1)';
yp_2.coefs(end-2:end) = (A3 \ by3)';

xp_2.coefs(4:end-3) = interp1([3 length(xp_2.coefs)-2],[xp_2.coefs(3) xp_2.coefs(end-2)],4:length(xp_2.coefs)-3);
yp_2.coefs(4:end-3) = interp1([3 length(yp_2.coefs)-2],[yp_2.coefs(3) yp_2.coefs(end-2)],4:length(yp_2.coefs)-3);

end

