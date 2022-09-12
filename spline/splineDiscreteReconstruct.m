function [breakPoints, splineValues, dsplineValues, ddsplineValues] = splineDiscreteReconstruct(knots, controlPoints)

breakPoints = unique(knots);
numCoefs = length(controlPoints);
cp = eye(numCoefs);

for i = 1:numCoefs
    sp = spmak(knots,cp(i,:));
    dsp = fnder(sp,1);
    ddsp = fnder(sp,2);
    for j = 1:length(breakPoints)
        if (i == 1)
            splineValues(1,j) = fnval(sp,breakPoints(j))*controlPoints(i);
            dsplineValues(1,j) = fnval(dsp,breakPoints(j))*controlPoints(i);
            ddsplineValues(1,j) = fnval(ddsp,breakPoints(j))*controlPoints(i);
        else
            splineValues(1,j) = splineValues(1,j) + fnval(sp,breakPoints(j))*controlPoints(i);
            dsplineValues(1,j) = dsplineValues(1,j) + fnval(dsp,breakPoints(j))*controlPoints(i);
            ddsplineValues(1,j) = ddsplineValues(1,j) + fnval(ddsp,breakPoints(j))*controlPoints(i);
        end
    end
end


end