function [value, isterminal, direction] = stickEvent(t,s,param,u,fCone,vec)

p = worldToObject(t,s,param,u);
value = checkSlip(p,fCone,vec);
isterminal = [1;1;1;1;1];
direction = [1;1;1;1;1];

end

