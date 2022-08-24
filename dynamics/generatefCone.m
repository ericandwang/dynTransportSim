function [fCone, vec] = generatefCone(param, maxNF)
% outputs 4 points corresponding to friction cone acceleration coordinates
% counterclockwise (positive n direction) in the object frame and normal
% vectors to the friciton cone facets

% Physical quantities
m_C = param(3);
I_C = param(4);
r_p1C = [param(5); param(6)];
r_p2C = [param(7); param(8)];
g = param(9);
mu = param(10);
J = [1 0 r_p1C(2); ...
    0 1 -r_p1C(1); ...
    1 0 r_p2C(2); ...
    0 1 -r_p2C(1)]; % Jacobian

% Query points in friction cone
P_query = zeros(4,4);
P_query(:,1) = [-mu*maxNF; maxNF; 0; 0];
P_query(:,2) = [mu*maxNF; maxNF; 0; 0];
P_query(:,3) = [0; 0; mu*maxNF; maxNF];
P_query(:,4) = [0; 0; -mu*maxNF; maxNF];

% Constructing generalized friction cone in object frame
W = J'*P_query;
A = [W(1,:)/m_C; W(2,:)/m_C; W(3,:)/I_C]; % acceleration space

% 4 point output (plane points)
fCone = A;

% vector output (plane normals)
vec = zeros(size(fCone));
for i = 1:4
    v1 = fCone(:,i);
    v2 = fCone(:,mod(i,4)+1);
    vec(:,i) = cross(v1,v2);
end

end

