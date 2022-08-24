function dc = dcGenTOPP_lin(r_GC, param, fCone, vec, dxp, ddxp, dyp, ddyp, ss, nPoints, th0)

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

for ii = 1:nPoints
    s = ss(ii);
    ds = dss(ii);
    dds = ddss(ii);
    th = ths(ii);
    dth = dths(ii);
    ddth = ddths(ii);
    ax_G = fnval(ddxp,s)*ds^2 + fnval(dxp,s)*dds;
    ay_G = fnval(ddyp,s)*ds^2 + fnval(dyp,s)*dds;
    p(:,ii) = [cos(th0(ii))-sin(th0(ii))*(th-th0(ii)) sin(th0(ii))+cos(th0(ii))*(th-th0(ii)) -R_GC*sin(th_GC); ...
              -(sin(th0(ii))+cos(th0(ii))*(th-th0(ii))) cos(th0(ii))-sin(th0(ii))*(th-th0(ii)) R_GC*cos(th_GC); ...
              0 0 1]*[ax_G; ay_G; ddth] + ...
              [g*(sin(th0(ii))+cos(th0(ii))*(th-th0(ii)))-R_GC*dth^2*cos(th_GC); ...
              g*(cos(th0(ii))-sin(th0(ii))*(th-th0(ii)))-R_GC*dth^2*sin(th_GC); 0];
end


% constructing linear plane constraints
dim = size(fCone,2);
for j = 1:nPoints
    for i = 1:dim
        c(i + dim*(j-1)) = vec(1,i)*(p(1,j)-fCone(1,i)) + vec(2,i)*(p(2,j)-fCone(2,i)) + ...
            vec(3,i)*(p(3,j)-fCone(3,i));
    end
end

jac = jacobian(c,states)';

% generating function for inequality gradient
dc = matlabFunction(jac,'Vars',{states});


end