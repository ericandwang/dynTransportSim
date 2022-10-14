function [f, gradf] = objTOPP(P,evalPoints)

% Objective function
% Integrating (summation) of inverse velocities
%f = sum(1./P(2:evalPoints-1,1));
f = sum(1./(P(1:evalPoints-1,1)+P(2:evalPoints,1)));

% Objective gradient
gradf = zeros(size(P));
gradf(1,1) = -1/(P(1,1) + P(2,1))^2;
gradf(evalPoints,1) = -1/(P(evalPoints-1,1) + P(evalPoints,1))^2;
gradf(2:evalPoints-1,1) = -1./(P(1:evalPoints-2,1)+P(2:evalPoints-1,1)).^2 - 1./(P(2:evalPoints-1,1)+P(3:evalPoints,1)).^2;
%gradf(2:evalPoints-1,1) = -1./(P(2:evalPoints-1,1).^2);

end

