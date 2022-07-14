function dists = checkSlip(p, fCone, vec)
% output dists gives orientation-dependent distances from p to each
% friction cone facet (positive outward facing normals)

nFacets = size(fCone,2);
dists = zeros(nFacets,1);
for i = 1:nFacets
    v = p - fCone(:,i);
    dists(i) = dot(v,vec(:,i))/norm(vec(:,i));
end
dists = [-p(2); dists];

end

