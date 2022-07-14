function [LfB, LgB, B] = barrierFuncGenJerk(r_GC, param, a, b)
% Generating barrier function Lie derivates and B(x_tilda) for jerk
% formulation

% parameters
g = param(9);
R_GC = norm(r_GC);
th_GC = angle(r_GC(1) + 1i*r_GC(2));

% symbolics
syms x dx y dy th dth ax_G ay_G ath_G

p = [cos(th) sin(th) -R_GC*sin(th_GC); ...
              -sin(th) cos(th) R_GC*cos(th_GC); ...
              0 0 1]*[ax_G; ay_G; ath_G] + ...
              [g*sin(th)-R_GC*dth^2*cos(th_GC); ...
              g*cos(th)-R_GC*dth^2*sin(th_GC); 0];

B_ = - p(1)^2/a - p(3)^2/b + p(2)^2;
s = [x; dx; y; dy; th; dth; ax_G; ay_G; ath_G];
A_dyn = [0 1 0 0 0 0 0 0 0; ...
         0 0 0 0 0 0 1 0 0; ...
         0 0 0 1 0 0 0 0 0; ...
         0 0 0 0 0 0 0 1 0; ...
         0 0 0 0 0 1 0 0 0; ...
         0 0 0 0 0 0 0 0 1; ...
         zeros(1,9); zeros(1,9); zeros(1,9)];
B_dyn = [zeros(6,3); ...
         eye(3)];
f_dyn = A_dyn*s;
g_dyn = B_dyn;

gradB = jacobian(B_,s);
LfB_ = gradB*f_dyn;
LgB_ = gradB*g_dyn;

% outputting functions from symbolics
B = matlabFunction(B_,'vars',s);
LfB = matlabFunction(LfB_,'vars',s);
LgB = matlabFunction(LgB_,'vars',s);

end