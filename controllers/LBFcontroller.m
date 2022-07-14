function [accel, comparison] = LBFcontroller(s0, param, aNominal, dT, u0, LfB, LgB, B, alpha, nu)
% LBF controller that outputs control input

% states
th = s0(5);
r_GC = [s0(7); s0(8)]; % object frame
    
% derived quantities
m = param(1) + param(3); % combined mass
m_C = param(3); % object mass
I = param(2) + param(4) + param(3)*norm(r_GC)^2; % combined inertia
g = param(9);
th_GC_ = angle(r_GC(1) + 1i*r_GC(2));

% optimization constraints
accel0 = (u0-[0;m*g;m_C*g*norm(r_GC)*cos(th_GC_+th)])./[m; m; I];
nonlincon_LBF_ = @(accel) nonlincon_LBF(dT,s0,accel0,accel,B,LfB,LgB,alpha,nu); % change for right nonlincon function

% optimization objective function
fun = @(accel) norm(accel-aNominal);
accel = fmincon(fun,accel0,[],[],[],[],[],[],nonlincon_LBF_);

% approximating change in barrier function
LfB_ = LfB(s0(1),s0(2),s0(3),s0(4),s0(5),s0(6),accel(1),accel(2),accel(3));
LgB_ = LgB(s0(1),s0(2),s0(3),s0(4),s0(5),s0(6),accel(1),accel(2),accel(3));
B_ = B(s0(1),s0(2),s0(3),s0(4),s0(5),s0(6),accel(1),accel(2),accel(3));
B_nom = B(s0(1),s0(2),s0(3),s0(4),s0(5),s0(6),accel0(1),accel0(2),accel0(3));
comparison = ((B_ - B_nom)/dT + LfB_ + LgB_*accel)*dT;


end

