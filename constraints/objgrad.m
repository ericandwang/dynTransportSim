function [f, gradf] = objgrad(x,x_g,nHorizon,costVector,finalCostVector)

x_goal = repmat([x_g;0;0;0],[nHorizon,1]);
costMatrix = diag([repmat(costVector,[nHorizon-1,1]); finalCostVector]);
f = (x-x_goal)'*costMatrix*(x-x_goal);
gradf = zeros(nHorizon*9,1);
for i = 1:nHorizon
    ind = 1+9*(i-1):9+9*(i-1);
    if i == nHorizon
        gradf(ind) = 2*finalCostVector.*(x(ind)) - 2*finalCostVector.*[x_g;0;0;0];
    else
        gradf(ind) = 2*costVector.*(x(ind)) - 2*costVector.*[x_g;0;0;0];
    end
end

end