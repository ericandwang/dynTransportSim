function [c, ceq] = nonlincon_LBF(dT, s0, accel0, accel, B, LfB, LgB, alpha, nu)

% defining Lie derivatives
LfB_ = LfB(s0(1),s0(2),s0(3),s0(4),s0(5),s0(6),accel(1),accel(2),accel(3));
LgB_ = LgB(s0(1),s0(2),s0(3),s0(4),s0(5),s0(6),accel(1),accel(2),accel(3));
B_ = B(s0(1),s0(2),s0(3),s0(4),s0(5),s0(6),accel(1),accel(2),accel(3));
B_nom = B(s0(1),s0(2),s0(3),s0(4),s0(5),s0(6),accel0(1),accel0(2),accel0(3));

% defining nonlinear inequality constraint
c = (B_ - B_nom)/dT + LfB_ + LgB_*accel + alpha*B_*sign(B_) - nu;
c = -c; % inverting so that <=

% defining nonlinear equality constraint
ceq = [];

end

