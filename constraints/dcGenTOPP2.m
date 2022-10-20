function dc = dcGenTOPP2(r_GC, param, fCone, vec, nPoints, accelLim)
% More efficient way of generating gradient UNUSED

syms ds dds th dth ddth dxs ddxs dys ddys
state = [ds; dds; th; dth; ddth];
spline = [dxs; ddxs; dys; ddys];

% Parameters
g = param(9);
R_GC = norm(r_GC);
th_GC = angle(r_GC(1) + 1i*r_GC(2));

ax_G = ddxs*ds^2 + dxs*dds;
ay_G = ddys*ds^2 + dys*dds;
p = [cos(th) sin(th) -R_GC*sin(th_GC); ...
     -sin(th) cos(th) R_GC*cos(th_GC); ...
     0 0 1]*[ax_G; ay_G; ddth] + ...
     [g*sin(th)-R_GC*dth^2*cos(th_GC); ...
     g*cos(th)-R_GC*dth^2*sin(th_GC); 0];

% constructing linear plane constraints
dim = size(fCone,2);
c1 = sym(zeros(dim,1));
for i = 1:dim
    c1(i) = vec(1,i)*(p(1)-fCone(1,i)) + vec(2,i)*(p(2)-fCone(2,i)) + ...
        vec(3,i)*(p(3)-fCone(3,i));
end
dc1 = jacobian(c1,state)';

% appending object frame y acceleration constraint
c2 = p(2) - accelLim;
dc2 = jacobian(c2,state)';

% generating intermediate functions
dc1_ = matlabFunction(dc1,'Vars',{state,spline});
dc2_ = matlabFunction(dc2,'Vars',{state,spline});

% concatenating all constraints
syms dss ddss ths dths ddths dxps ddxps dyps ddyps
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

for i = 1:nPoints
    stateI = [dss(i); ddss(i); ths(i); dths(i); ddths(i)];
    splineI = [dxps(i); ddxps(i); dyps(i); ddyps(i)];
    dc1s{i} = dc1_(stateI,splineI);
    dc2s{i} = dc2_(stateI,splineI);
end

% create block diagonal matrices
dc1 = blkdiag(dc1s{:});
dc2 = blkdiag(dc2s{:});
dcs = [dc1, dc2];

% Rearrange rows so that reflects order of states
dc_ = dcs;
for i = 1:5
    dc_(1+nPoints*(i-1):nPoints*i,:) = dcs(i:5:end,:);
end

dc = matlabFunction(dc_,'Vars',{states,splines});


end