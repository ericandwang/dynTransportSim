function [forces] = individualForces(vec, plist)

forces = zeros(4,length(plist)); % each column: [fxl; fyl; fxr; fyr]
nl = vec(:,1)/norm(vec(:,1));
nr = vec(:,3)/norm(vec(:,3));

for i = 1:length(plist)
    p = plist(:,i);
    pl = [p(1:2); -(nl(1)*p(1)+nl(2)*p(2))/nl(3)];
    pr = [p(1:2); -(nr(1)*p(1)+nr(2)*p(2))/nr(3)];
    A = [pl(2:3), pr(2:3)];
    scale = A \ p(2:3);
    forces(:,i) = [p(1:2)*scale(1); p(1:2)*scale(2)];
end

end

