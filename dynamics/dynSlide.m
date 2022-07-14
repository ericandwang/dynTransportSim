function ds = dynSlide(t,s,param,u)
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
dr_GC = [s(9); s(10)]; % object frame

% parameters
m_G = param(1);
I_G = param(2);
m_C = param(3);
I_C = param(4);
%r_p1C = [param(5); param(6)]; % object frame
%r_p2C = [param(7); param(8)]; % object frame
g = param(9);
mu = param(10);

% control inputs
Fx = u(1);
Fy = u(2);
Fth = u(3);

% derived quantities
r = r_GC(1); dr = dr_GC(1);
h = r_GC(2); dh = dr_GC(2);

% system of equations for accelerations
% accelerations = [ddx_G; ddy_G; ddth; ddr; N; M]
A = [cos(th) sin(th) -h 1 mu*sign(dr) 0; ...
     -sin(th) cos(th) r 0 -1 0; ...
     0 0 I_G 0 r+mu*sign(dr)*h 1; ...
     0 0 I_C 0 0 -1; ...
     m_G 0 0 0 -mu*sign(dr)*cos(th)-sin(th) 0; ...
     0 m_G 0 0 -mu*sign(dr)*sin(th)+cos(th) 0];
b = [r*dth^2-m_C*g*sin(th); ...
     -2*dr*dth+h*dth^2-m_C*g*cos(th); ...
     Fth; ...
     0; ...
     Fx; ...
     Fy-m_G*g];
accelerations = A\b;

% dynamics
ddx_G = accelerations(1);
ddy_G = accelerations(2);
ddth = accelerations(3);
ddr = accelerations(4);
ddr_GC = [ddr; 0];

% output change of state vector
ds = [dx_G; ...
    ddx_G; ...
    dy_G; ...
    ddy_G; ...
    dth; ...
    ddth; ...
    dr_GC; ...
    ddr_GC];

end
