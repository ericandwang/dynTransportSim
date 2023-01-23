function [f, gradf] = objTOPPLambda(P)

% Objective function
% Max lambda
f = -P(end);

% Objective gradient
gradf = zeros(size(P));
gradf(end) = -1;

end

