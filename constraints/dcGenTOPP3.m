function fun = dcGenTOPP3(r_GC, param, fCone, vec, nPoints, accelLim)
% gen function

syms dss ddss ths dths ddths
dss = sym('dss',[nPoints,1]);
ddss = sym('ddss',[nPoints,1]);
ths = sym('ths',[nPoints,1]);
dths = sym('dths',[nPoints,1]);
ddths = sym('ddths',[nPoints,1]);
dxps = sym('dxps',[nPoints,1]);
ddxps = sym('ddxps',[nPoints,1]);
dyps = sym('dyps',[nPoints,1]);
ddyps = sym('ddyps',[nPoints,1]);
states = [dss; ddss; ths; dths; ddths];
splines = [dxps; ddxps; dyps; ddyps];

% Parameters
g = param(9);
R_GC = norm(r_GC);
th_GC = angle(r_GC(1) + 1i*r_GC(2));

for ii = 1:nPoints
    th = ths(ii);
    dth = dths(ii);
    ddth = ddths(ii);
    ax_G = ddxps(ii)*dss(ii)^2 + dxps(ii)*ddss(ii);
    ay_G = ddyps(ii)*dss(ii)^2 + dyps(ii)*ddss(ii);
    p(:,ii) = [cos(th) sin(th) -R_GC*sin(th_GC); ...
              -sin(th) cos(th) R_GC*cos(th_GC); ...
              0 0 1]*[ax_G; ay_G; ddth] + ...
              [g*sin(th)-R_GC*dth^2*cos(th_GC); ...
              g*cos(th)-R_GC*dth^2*sin(th_GC); 0];
end

% constructing linear plane constraints  (knot points and collocation points)
nPoints = size(p,2);
dim = size(fCone,2);
for j = 1:nPoints
    for i = 1:dim
        c(i + dim*(j-1),1) = vec(1,i)*(p(1,j)-fCone(1,i)) + vec(2,i)*(p(2,j)-fCone(2,i)) + ...
            vec(3,i)*(p(3,j)-fCone(3,i));
    end
end

% appending object frame y acceleration
cappend = p(2,:)' - accelLim;
c = [c; cappend];

jac = jacobian(c,states)';

% generating function for inequality gradient
if (nargout == 0)
    matlabFunction(jac,'Vars',{states,splines},'File','constraints/gen/dcFunGen.m');
else
    fun = matlabFunction(jac,'Vars',{states,splines});
end


end