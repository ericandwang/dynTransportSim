function [c, ceq] = nonlinconPlaneCBF(nHorizon, x, s, param, Bp, Bp_initial, gamma)

% initial states
r_GC = [s(7); s(8)]; % object frame

% parameters
g = param(9);

% defining vector of barrier functions
Bvec_ = zeros(4*(nHorizon+1),1);
gammavec = repmat(gamma,nHorizon,1);
Bvec_(1:4,1) = Bp_initial;
for i = 1:nHorizon
    s_u = x(1+9*(i-1):9+9*(i-1));
    Bvec_(1+4*i:4+4*i,1) = Bp(s_u);
end

% constructing barrier function nonlinear inequalities
c = Bvec_(1:end-4).*(1-gammavec) - Bvec_(5:end);

% no nonlinear equality constraints
ceq = [];

end

