function dceq = dceqGenPATH(numCoefs, velEndpoints, s0, dx_des)

% spline construction
syms coefs
coefs = sym('coefs',[numCoefs,1]);
numCoefs = length(coefs);
xcoefs = coefs(1:numCoefs/2);
ycoefs = coefs(numCoefs/2+1:end);

vxEndpoints = velEndpoints*xcoefs; % actual spline endpoint derivatives
vyEndpoints = velEndpoints*ycoefs;
vxDes = [s0(2); dx_des(1)]; % desired velocities
vyDes = [s0(4); dx_des(2)];

% initial
if (vxDes(1) == 0 && vyDes(1) == 0)
    ceq1(1,1) = vxEndpoints(1);
    ceq1(2,1) = vyEndpoints(1);
else
    if (vxDes(1) == 0)
        ceq1 = vxEndpoints(1);
    elseif (vyDes(1) == 0)
        ceq1 = vyEndpoints(1);
    else
        ceq1 = vxEndpoints(1)/vyEndpoints(1) - vxDes(1)/vyDes(1);
    end
end
% final
if (vxDes(2) == 0 && vyDes(2) == 0)
    ceq2(1,1) = vxEndpoints(2);
    ceq2(2,1) = vyEndpoints(2);
else
    if (vxDes(2) == 0)
        ceq2 = vxEndpoints(2);
    elseif (vyDes(2) == 0)
        ceq2 = vyEndpoints(2);
    else
        ceq2 = vxEndpoints(2)/vyEndpoints(2) - vxDes(2)/vyDes(2);
    end
end
ceq = [ceq1; ceq2];


jac = jacobian(ceq,coefs)';

% generating function for inequality gradient
dceq = matlabFunction(jac,'Vars',{coefs});

end

