function [accel, comparison] = LBFcontrollerJerk(s0, param, aNominal, dT, u0, LfB, LgB, B, alpha, LfB_reg, LgB_reg, B_reg)
% LBF controller that outputs control input based on jerk integration

    
    % states
    th = s0(5);
    r_GC = [s0(7); s0(8)]; % object frame
    
    % derived quantities
    m = param(1) + param(3); % combined mass
    m_C = param(3); % object mass
    I = param(2) + param(4) + param(3)*norm(r_GC)^2; % combined inertia
    g = param(9);
    th_GC_ = angle(r_GC(1) + 1i*r_GC(2));
    
    % evaluating symbolic functions
    accel0 = (u0-[0;m*g;m_C*g*norm(r_GC)*cos(th_GC_+th)])./[m; m; I];
    LfB_ = LfB(s0(1),s0(2),s0(3),s0(4),s0(5),s0(6),accel0(1),accel0(2),accel0(3));
    LgB_ = LgB(s0(1),s0(2),s0(3),s0(4),s0(5),s0(6),accel0(1),accel0(2),accel0(3));
    B_ = B(s0(1),s0(2),s0(3),s0(4),s0(5),s0(6),accel0(1),accel0(2),accel0(3));
    
    H = eye(3).*dT^2; % using for quadprog
    f = (accel0-aNominal)*dT;
    C = eye(3).*dT; % using for lsqnonlin
    d = (aNominal-accel0);
    A = -LgB_;
    b = LfB_ + alpha*B_*sign(B_);
    
    %jerk0 = quadprog(H,f); % using for quadprog
    %jerk = quadprog(H,f,A,b,[],[],jerk0);
    %jerkTest = quadprog(H,f,A,b);
    opts1=  optimset('display','off');
    jerkTest = lsqlin(C,d,A,b,[],[],[],[],[],opts1);
    
    
    % searching along jerk and finding correct dT in order to guarantee CBF
    % condition met
    %n = 1000000;
    %subd = 100000;
    n = 1000;
    subd = dT/n;
    Bplot = zeros(fix(n),1);
    Bdotplot_reg = zeros(fix(n),1);
    for i = 1:n
        accelNew = accel0 + jerkTest*i*subd;
        Bplot(i) = B(s0(1),s0(2),s0(3),s0(4),s0(5),s0(6),accelNew(1),accelNew(2),accelNew(3))-B_;
        LfB_reg_ = LfB_reg(s0(1),s0(2),s0(3),s0(4),s0(5),s0(6),accelNew(1),accelNew(2),accelNew(3));
        LgB_reg_ = LgB_reg(s0(1),s0(2),s0(3),s0(4),s0(5),s0(6),accelNew(1),accelNew(2),accelNew(3));
        B_reg_ = B_reg(s0(1),s0(2),s0(3),s0(4),s0(5),s0(6),accelNew(1),accelNew(2),accelNew(3));
        Bdotplot_reg(i) = LfB_reg_ + LgB_reg_*accelNew;
    end
    Bd = Bplot + Bdotplot_reg*dT; % added Bdotplot_reg*dT to predict Bdot contribution from subsequent state evolution
    % close all
    % figure, plot(Bplot)
    % yyaxis right
    % hold on, plot(Bdotplot_reg);
    % figure, plot(Bd)
    % hold on, plot(B_ + Bd)
    [~,ind] = max(Bd);
    comparison = Bplot(ind) + Bdotplot_reg(ind)*dT;
    jerk = jerkTest;
    dT = ind*subd;
    
    accel = accel0 + jerk.*dT;


end

