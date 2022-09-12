function [coefs, dcoefs, ddcoefs] = splineBasisCoefs(knots, order)

numBasis = length(knots) - order;
cp = eye(numBasis);
coefs = cell(numBasis,1);
ddImpulses = zeros(numBasis,numBasis-2);

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
end

end

