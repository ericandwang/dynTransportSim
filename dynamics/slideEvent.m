function [value, isterminal, direction] = slideEvent(t,s,param,u,fCone,vec)

% parameters
handW = param(11);
handH = param(12);
objectW = param(13);
objectH = param(14);

% states
r_GC = [s(7); s(8)]; % object frame
dr_GC = [s(9); s(10)]; % object frame
r = r_GC(1); dr = dr_GC(1);

value = [r+objectW/2-handW/2; r-objectW/2+handW/2; dr];
isterminal = [1;1;1];
direction = [1;-1;0];

end