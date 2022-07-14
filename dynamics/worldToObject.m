function p = worldToObject(t,s,param,u)

% states
th = s(5);
dth = s(6);
r_GC = [s(7); s(8)]; % object frame

% dynamics
ds = dynStick(t,s,param,u);

% rotation matrix
%rot = @(th) [cos(th) -sin(th); sin(th) cos(th)];

% converting acceleration at G from world frame to object frame (_)
% world accelerations
ax_G = ds(2);
ay_G = ds(4);
ath_G = ds(6);

% parameters
g = param(9);

% converting to object acceleration p
R_GC = norm(r_GC);
th_GC = angle(r_GC(1) + 1i*r_GC(2));
p = [cos(th) sin(th) -R_GC*sin(th_GC); ...
 -sin(th) cos(th) R_GC*cos(th_GC); ...
 0 0 1]*[ax_G; ay_G; ath_G] + ...
    [g*sin(th)-R_GC*dth^2*cos(th_GC); ...
     g*cos(th)-R_GC*dth^2*sin(th_GC); 0];


% Previous code
% tempVec = rot(-th)*[ax_G;ay_G];
% ax_G_ = tempVec(1); ay_G_ = tempVec(2); ath_G_ = ath_G;
% centrifugal = -dth^2*[r_GC;0];
% tangential = cross([0;0;ath_G_],[r_GC;0]);
% a_ = [ax_G_;ay_G_;ath_G_] + centrifugal + tangential;
% 
% % parameters
% m_G = param(1);
% I_G = param(2);
% m_C = param(3);
% I_C = param(4);
% g = param(9);
% g_ = [rot(-th)*[0;g]; 0];
% 
% % control inputs
% %Fx = u(1);
% %Fy = u(2);
% %Fth = u(3);
% 
% % derived quantities
% %m = m_G + m_C; % combined mass
% %I = I_G + I_C + m_C*norm(r_GC)^2; % combined inertia
% 
% p = a_ + g_; % acceleration point query in generalized friction cone

end

