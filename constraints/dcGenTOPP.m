function dc = dcGenTOPP(r_GC, param, fCone, vec, xp, dxp, ddxp, yp, dyp, ddyp, ss, nPoints, accelLim)

syms dss ddss ths dths ddths
dss = sym('dss',[nPoints,1]);
ddss = sym('ddss',[nPoints,1]);
ths = sym('ths',[nPoints,1]);
dths = sym('dths',[nPoints,1]);
ddths = sym('ddths',[nPoints,1]);
states = [dss; ddss; ths; dths; ddths];

% DC = cell(nPoints,1);
% % inequality gradients
% for ii = 1:nPoints
%     jac = jacob(ss(ii));
%     DC{ii} = jac([dss(ii); ddss(ii); ths(ii); dths(ii); ddths(ii)]);
% end
% dcsym = blkdiag(DC{:});
% dcreorder = dcsym;
% for ii = 1:nPoints
%     dcreorder(ii,:) = dcsym(1+5*(ii-1),:);
%     dcreorder(ii + nPoints,:) = dcsym(2+5*(ii-1),:);
%     dcreorder(ii + 2*nPoints,:) = dcsym(3+5*(ii-1),:);
%     dcreorder(ii + 3*nPoints,:) = dcsym(4+5*(ii-1),:);
%     dcreorder(ii + 4*nPoints,:) = dcsym(5+5*(ii-1),:);
% end
% 
% % generating function for inequality gradient
% dc = matlabFunction(dcreorder,'Vars',{states});

% Parameters
g = param(9);
R_GC = norm(r_GC);
th_GC = angle(r_GC(1) + 1i*r_GC(2));

xp_ = fnval(xp,ss);
dxp_ = fnval(dxp,ss);
ddxp_ = fnval(ddxp,ss);
yp_ = fnval(yp,ss);
dyp_ = fnval(dyp,ss);
ddyp_ = fnval(ddyp,ss);
for ii = 1:nPoints
    ds = dss(ii);
    dds = ddss(ii);
    th = ths(ii);
    dth = dths(ii);
    ddth = ddths(ii);
    ax_G = ddxp_(ii)*ds^2 + dxp_(ii)*dds;
    ay_G = ddyp_(ii)*ds^2 + dyp_(ii)*dds;
    p(:,ii) = [cos(th) sin(th) -R_GC*sin(th_GC); ...
              -sin(th) cos(th) R_GC*cos(th_GC); ...
              0 0 1]*[ax_G; ay_G; ddth] + ...
              [g*sin(th)-R_GC*dth^2*cos(th_GC); ...
              g*cos(th)-R_GC*dth^2*sin(th_GC); 0];
end

% appending collocation points CCC unused
% h = ss(2)-ss(1);
% for ii = 1:nPoints-1
%     dds = ddss(ii);
%     ds = sqrt(h*dds + dss(ii)^2);
%     dx = -3/(2*h)*(xp_(ii)-xp_(ii+1)) - 1/4*(dxp_(ii)+dxp_(ii+1));
%     dy = -3/(2*h)*(yp_(ii)-yp_(ii+1)) - 1/4*(dyp_(ii)+dyp_(ii+1));
%     ddx = (ddxp_(ii)+ddxp_(ii+1))/2;
%     ddy = (ddyp_(ii)+ddyp_(ii+1))/2;
%     dt = h/(dss(ii)+ds);
%     ddth = ddths(ii);
%     dth = dths(ii) + ddth*dt;
%     th = ths(ii) + dth*dt + 1/2*ddth*dt^2;
%     ax_G = ddx*ds^2 + dx*dds;
%     ay_G = ddy*ds^2 + dy*dds;
%     pappend(:,ii) = [cos(th) sin(th) -R_GC*sin(th_GC); ...
%               -sin(th) cos(th) R_GC*cos(th_GC); ...
%               0 0 1]*[ax_G; ay_G; ddth] + ...
%               [g*sin(th)-R_GC*dth^2*cos(th_GC); ...
%               g*cos(th)-R_GC*dth^2*sin(th_GC); 0];
% end

%p = [p, pappend];


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
dc = matlabFunction(jac,'Vars',{states});


end