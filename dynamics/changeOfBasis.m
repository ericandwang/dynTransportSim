function [basisVectors, constraintSlopes] = changeOfBasis(fCone)

th = atan((fCone(3,2) - fCone(3,1))/(fCone(1,2) - fCone(1,1)));
basisVectors = [cos(th) 0 -sin(th); ...
                0 1 0; ...
                sin(th) 0  cos(th)];
fCone_ = basisVectors \ fCone;
constraintSlopes(1) = abs(fCone(2,1))/abs(fCone(1,1));
constraintSlopes(2) = abs(fCone(2,1))/abs(fCone(3,1));
% p_ = basisVectors \ p; converting cartesian coordinate to rotated basis
                       % coordinates

end