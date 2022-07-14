function accel = impedenceController(s0, param, sNominal, uNominal, gains)
% Simple impedence controller that has a feedforward term based on nominal
% control and a feedback term based on hand-tuned gains

% states
r_GC = [s0(7); s0(8)]; % object frame

% derived quantities
m = param(1) + param(3); % combined mass
I = param(2) + param(4) + param(3)*norm(r_GC)^2; % combined inertia

x = s0(1);
y = s0(3);
th = s0(5);
x_nom = sNominal(1);
y_nom = sNominal(3);
th_nom = sNominal(5);

ff = uNominal./[m; m; I];
fb = gains.*[x_nom-x; y_nom-y; th_nom-th];

accel = ff + fb;

end

