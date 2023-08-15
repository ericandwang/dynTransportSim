function initialVel = sVelStatic(xp, vx, sBounds, evalPoints)
% inputs
% xp (x position spline)
% yp (y position spline)
% vx (2x1 array [vx0, vxf])

% outputs
% initialVel (s values for 1...evalPoints spline evaluation points given
% constant acceleration)

evalPoints = 10000;

initialVel = ones(evalPoints,1);
s = linspace(sBounds(1),sBounds(2),evalPoints);
a = (vx(2)^2 - vx(1)^2)/(2*(fnval(xp,sBounds(2))-fnval(xp,sBounds(1))));
for i = 2:evalPoints-1
    v = sign(vx(1)+vx(2))*sqrt(vx(1)^2 + 2*a*(fnval(xp,s(i))-fnval(xp,sBounds(1))));
    initialVel(i) = v/fnval(fnder(xp,1),s(i));
end
if (vx(1) ~= 0)
    initialVel(1) = vx(1)/fnval(fnder(xp,1),s(1));
else
    initialVel(1) = initialVel(2);
end
if (vx(2) ~= 0)
    initialVel(end) = vx(2)/fnval(fnder(xp,1),s(end));
else
    initialVel(end) = initialVel(end-1);
end

% debugging 8/15
% v = initialVel.*fnval(fnder(xp,1),s)';
% 
% dss0 = ones(evalPoints,1).*initialVel;
% ddss0 = zeros(evalPoints,1);
% for ii = 1:length(ddss0)-1
%     ddss0(ii) = (dss0(ii+1)^2 - dss0(ii)^2)/(2*(s(ii+1)-s(ii)));
% end
% for i = 1:evalPoints
%     ax(i) = fnval(fnder(xp,2),s(i))*dss0(i)^2 + fnval(fnder(xp,1),s(i))*ddss0(i);
% end

end

