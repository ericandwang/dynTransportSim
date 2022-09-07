function coefs = splineBasisCoefs(knots, order)

numBasis = length(knots) - order;
cp = eye(numBasis);
coefs = cell(numBasis,1);

for i = 1:numBasis
    bs = spmak(knots,cp(i,:));
    sp = fn2fm(bs,'pp');
    coefs{i} = sp.coefs;
end

end

