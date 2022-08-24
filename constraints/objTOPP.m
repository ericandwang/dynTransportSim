function [f, gradf] = objTOPP(P,evalPoints)

% Objective function
% Integrating (summation) of inverse velocities
f = sum(1./P(2:evalPoints-1,1));

% Objective gradient
gradf = zeros(size(P));
gradf(2:evalPoints-1,1) = -1./(P(2:evalPoints-1,1).^2);

end
