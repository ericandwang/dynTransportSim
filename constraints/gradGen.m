function DC = gradGen(r_GC, param, fCone, vec)

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

% calculating nonlinear constraint form
dim = size(fCone,2);
for i = 1:dim
    c(i,1) = vec(1,i)*(p(1)-fCone(1,i)) + vec(2,i)*(p(2)-fCone(2,i)) + ...
            vec(3,i)*(p(3)-fCone(3,i));
end

% calculating jacobian
s = [x; dx; y; dy; th; dth; ax_G; ay_G; ath_G];
jac = jacobian(c,s)';

% generating function for jacobian
DC = matlabFunction(jac,'vars',s);

end