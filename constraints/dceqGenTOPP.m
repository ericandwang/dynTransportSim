function dceq = dceqGenTOPP(ss, nPoints)

%syms dss ddss ths dths ddths
dss = sym('dss',[nPoints,1]);
ddss = sym('ddss',[nPoints,1]);
ths = sym('ths',[nPoints,1]);
dths = sym('dths',[nPoints,1]);
ddths = sym('ddths',[nPoints,1]);
states = [dss; ddss; ths; dths; ddths];

% path acceleration constraints
ceq1 = dss(2:end).^2 - dss(1:end-1).^2 - ...
    2.*(ss(2:end)-ss(1:end-1)).*ddss(1:end-1);

% orientation kinematic constraints
dt = 2.*(ss(2:end)-ss(1:end-1))./(dss(2:end)+dss(1:end-1));
ceq2 = dths(2:end)-dths(1:end-1) - ddths(1:end-1).*dt;
ceq3 = (ths(1:end-1)-ths(2:end)) + dths(1:end-1).*dt + ...
    1/2.*ddths(1:end-1).*dt.^2;

ceq = [ceq1; ceq2; ceq3];
jac = jacobian(ceq,states)';

% generating function for equality gradient
dceq = matlabFunction(jac,'Vars',{states});

end