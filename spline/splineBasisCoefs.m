function [coefs, dcoefs, ddcoefs, posEndpoints, velEndpoints, accelEndpoints] = splineBasisCoefs(knots, order)

numBasis = length(knots) - order;
cp = eye(numBasis);
coefs = cell(numBasis,1);
%ddImpulses = zeros(numBasis,numBasis-2);
posEndpoints = size(2,numBasis);
velEndpoints = size(2,numBasis);
accelEndpoints = size(2,numBasis);

for i = 1:numBasis
    bs = spmak(knots,cp(i,:));
    dbs = fnder(bs,1);
    ddbs = fnder(bs,2);
    sp = fn2fm(bs,'pp');
    dsp = fn2fm(dbs,'pp');
    ddsp = fn2fm(ddbs,'pp');
    coefs{i} = sp.coefs;
    dcoefs{i} = dsp.coefs;
    ddcoefs{i} = ddsp.coefs;
    posEndpoints(1,i) = fnval(bs,knots(1));
    posEndpoints(2,i) = fnval(bs,knots(end));
    velEndpoints(1,i) = fnval(dbs,knots(1));
    velEndpoints(2,i) = fnval(dbs,knots(end));
    accelEndpoints(1,i) = fnval(ddbs,knots(1));
    accelEndpoints(2,i) = fnval(ddbs,knots(end));
end

end

