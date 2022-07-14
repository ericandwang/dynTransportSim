function ds = dynStick(t,s,param,u)
% Note: unless specified in the object frame, all quantities are in world
% frame

% states
x_G = s(1);
dx_G = s(2);
y_G = s(3);
dy_G = s(4);
th = s(5);
dth = s(6);
r_GC = [s(7); s(8)]; % object frame

% parameters
m_G = param(1);
I_G = param(2);
m_C = param(3);
I_C = param(4);
%r_p1C = [param(5); param(6)]; % object frame
%r_p2C = [param(7); param(8)]; % object frame
g = param(9);

% control inputs
Fx = u(1);
Fy = u(2);
Fth = u(3);

% derived quantities
m = m_G + m_C; % combined mass
I = I_G + I_C + m_C*norm(r_GC)^2; % combined inertia
th_GC_ = atan(r_GC(2)/r_GC(1)); % angle between r_GC and object frame x axis

% dynamics
ddx_G = Fx/m;
ddy_G = (Fy - m*g)/m;
ddth = (Fth - m_C*g*norm(r_GC)*cos(th_GC_+th))/I;
dr_GC = [0; 0];

% output change of state vector
ds = [dx_G; ...
    ddx_G; ...
    dy_G; ...
    ddy_G; ...
    dth; ...
    ddth; ...
    dr_GC];

end
