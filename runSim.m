close all
clear all
folder = '../dynTransportSim';
addpath(genpath(folder))

%% Parameter Setup
% hand and object dimensions
handW = 0.26*2*3;
handH = 0.01;
objectW = 0.10;
objectH = 0.10*1;
% parameters
m_G = 1.2; % kg
I_G = 0.1; % kg*m^2
m_C = 0.4; % kg
I_C = 0.05; % kg*m^2
r_p1C = [objectW/2; objectH/2]; % m
r_p2C = [-objectW/2; objectH/2]; % m
g = 9.81;
mu = 0.7;%0.1; %0.3 %0.5
param = [m_G; I_G; m_C; I_C; r_p1C; r_p2C; g; mu; handW; handH; objectW; objectH];
% controller actuation limits (bidirectional)
u_max = [30; 30; 20];
% desired position
g_des = [10; 10];
% tolerance for distance from plane constraints
tol = 0.05;%0.01 (0.05 = conservative);
% Parameters for control barrier function cone
a = (2.4/25)^2; % most conservative (2.5/25)^2
b = (8.9/25)^2; % most conservative (9/25)^2

%% Generate Friction Cone (object frame)
maxNormalForce = 10; % max normal force [N]
[fCone, vec] = generatefCone(param,maxNormalForce);

% predicted fCone and vec
fCone_predict = fCone;
%fCone_predict(1,1:4) = fCone_predict(1,1:4).*2.56./2.5; % mu = 0.1
%fCone_predict(1,1:4) = fCone_predict(1,1:4).*12.58./12.5; % mu = 0.5
fCone_predict(1,1:4) = fCone_predict(1,1:4).*17.6./17.5; % mu = 0.7
vec_predict = zeros(size(fCone_predict));
for i = 1:4
    v1 = fCone_predict(:,i);
    v2 = fCone_predict(:,mod(i,4)+1);
    vec_predict(:,i) = cross(v1,v2);
end

trippedSlip = false;

%fCone_predict = fCone;
%vec_predict = vec;


%% Simulation Setup
% simulation time steps
controllerBandwidth = 100; % Hz
deltaT = 1/controllerBandwidth; % sec
tDuration = 10; % sec
timeSteps = tDuration/deltaT;
% controller options
stabilizingControl = [3 6 7 8 9];
saveController = 0;
fileName = 'trajOpt3-29.mat';
saveVideo = 0;
specialName = '';
controlMode = 10;
mode = 0; % -1 = exit (fail) 0 = dynamic grasping 1 = sliding
if (any(stabilizingControl == controlMode))
    load(fileName)
end
% initial conditions
if (any(stabilizingControl == controlMode))
    %disturbance = zeros(8,1);
    disturbance = [-0.2; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0]; %-0.2 x
    s0 = s0Save + disturbance;
    r_GC = [s0(7); s0(8)]; % object frame
    th_GC_ = atan(r_GC(2)/r_GC(1)); % angle between r_GC and object frame x axis
    u0 = [0; (m_G+m_C)*g; m_C*g*norm(r_GC)*cos(th_GC_+s0(5))]; %CCC fix by saving u0
else
    s0 = [-5; ... % x_G -5
        0; ...   % dx_G
        0; ...   % y_G 5
        0; ...   % dy_G
        0; ...   % th
        0; ...   % dth
        0.07; ...   % r_GC 0.07
        objectH/2+handH/2; ...
        0; ... %dr_GC
        0; ...
        0; ... %th_C
        0]; %dth_C
    r_GC = [s0(7); s0(8)]; % object frame
    th_GC_ = atan(r_GC(2)/r_GC(1)); % angle between r_GC and object frame x axis
    u0 = [0; (m_G+m_C)*g; m_C*g*norm(r_GC)*cos(th_GC_+s0(5))];
end

%% Calculate Constraint Jacobian
r_GC = [s0(7);s0(8)];
jacob = gradGen(r_GC, param, fCone, vec);

%% Calculate LBF functions
[LfB_J, LgB_J, B_J] = barrierFuncGenJerk(r_GC, param, a, b);
[LfB, LgB, B] = barrierFuncGen(r_GC, param, a, b);
[LfB1, LgB1, B1] = barrierFuncGenPlane(r_GC, param, 1, fCone_predict, vec_predict, tol);
[LfB2, LgB2, B2] = barrierFuncGenPlane(r_GC, param, 2, fCone_predict, vec_predict, tol);
[LfB3, LgB3, B3] = barrierFuncGenPlane(r_GC, param, 3, fCone_predict, vec_predict, tol);
[LfB4, LgB4, B4] = barrierFuncGenPlane(r_GC, param, 4, fCone_predict, vec_predict, tol);
Bp = @(s_u) [B1(s_u); B2(s_u); B3(s_u); B4(s_u)];

%% Run Simulation
tTraj = [];
sTraj = [];
uTraj = [];
mTraj = [];
u = u0;
comparison = [];
currentTime = 0;
if (saveController)
    uSave = [];
    s0Save = s0;
    tSave = [];
    sSave = [];
end
for i = 1:timeSteps
    % control input
    if (controlMode == 1) % 1. gravity compensation
        % states
        th = s0(5);
        r_GC = [s0(7); s0(8)]; % object frame

        % derived quantities
        m = m_G + m_C; % combined mass
        I = I_G + I_C + m_C*norm(r_GC)^2; % combined inertia
        th_GC_ = atan(r_GC(2)/r_GC(1)); % angle between r_GC and object frame x axis
        if (currentTime < 2.23)
            u = [0; (m_G+m_C)*g; m_C*g*norm(r_GC)*cos(th_GC_+th) - (th-pi/8)/10];
        else
            u = [0; (m_G+m_C)*g; m_C*g*norm(r_GC)*cos(th_GC_+th) - 2*3*(th+pi/20)/10];
        end
    elseif (controlMode == 2) % 2. random sampling
        outFrictionCone = true;
        while (outFrictionCone)
            u = 2*rand(3,1).*u_max-u_max;
            [value, ~, ~] = stickEvent((i-1)*deltaT,s0,param,u,fCone,vec);
            if all(value < -0.5)
                outFrictionCone = false;
            end
        end
    elseif (controlMode == 3) % 3. used saved open loop policy
        try
            u = uSave(i,:)';
        catch
            u = zeros(3,1);
        end
    elseif (controlMode == 4) % 4. use single step controller
        % states
        th = s0(5);
        r_GC = [s0(7); s0(8)]; % object frame

        % derived quantities
        m = m_G + m_C; % combined mass
        I = I_G + I_C + m_C*norm(r_GC)^2; % combined inertia
        th_GC_ = atan(r_GC(2)/r_GC(1)); % angle between r_GC and object frame x axis
        
        pGains = [2; 2; 2];
        dGains = [0.1; 0.1; 0.1]*5;
        a_des = pGains.*[-s0(1); -s0(3); -s0(5)] + ...
            dGains.*[-s0(2); -s0(4); -s0(6)]; % simple desired acceleration > positional error
        accel = controller(deltaT, 1, s0, param, fCone, vec, tol, jacob, a_des);
        u = [m*accel(1); m*accel(2)+m*g; I*accel(3)+m_C*g*norm(r_GC)*cos(th_GC_+th)];
    elseif (controlMode == 5) % 5. use receding horizon controller
        % states
        th = s0(5);
        r_GC = [s0(7); s0(8)]; % object frame

        % derived quantities
        m = m_G + m_C; % combined mass
        I = I_G + I_C + m_C*norm(r_GC)^2; % combined inertia
        th_GC_ = atan(r_GC(2)/r_GC(1)); % angle between r_GC and object frame x axis
        
        a_des = [0; 0; 0]; % desired position is the origin
        nHorizon = 10; % number of steps
        dt = deltaT*10; % time delta between each step
        x_g = [0;0;0;0;0;0];
        accel = controller(dt, nHorizon, s0, param, fCone, vec, tol, jacob, a_des, x_g);
        u = [m*accel(1); m*accel(2)+m*g; I*accel(3)+m_C*g*norm(r_GC)*cos(th_GC_+th)];
    elseif (controlMode == 6) % 6. use impedence controller to stabilize trajectory
        % states
        r_GC = [s0(7); s0(8)]; % object frame

        % derived quantities
        m = m_G + m_C; % combined mass
        I = I_G + I_C + m_C*norm(r_GC)^2; % combined inertia

        [~, ind] = min(abs(deltaT*i-tSave));
        sNominal = sSave(ind,:)';
        uNominal = uSave(i,:)';
        gains = [1; 1; 10]; % x y th
        accel = impedenceController(s0,param,sNominal,uNominal,gains);
        u = [m*accel(1); m*accel(2); I*accel(3)];
    elseif (controlMode == 7) % 7. LBF controller (jerk formulation)
        % states
        th = s0(5);
        r_GC = [s0(7); s0(8)]; % object frame

        % derived quantities
        m = m_G + m_C; % combined mass
        I = I_G + I_C + m_C*norm(r_GC)^2; % combined inertia
        th_GC_ = atan(r_GC(2)/r_GC(1)); % angle between r_GC and object frame x axis

        [~, ind] = min(abs(deltaT*i-tSave));
        sNominal = sSave(ind,:)';
        uNominal = uSave(i,:)';
        gains = [1; 1; 10]; % x y th
        alpha = -1;%10e4;
        accel_impedence = impedenceController(s0,param,sNominal,uNominal,gains);
        accel_intermediate = accel_impedence-[0;g;0];
        accel_intermediate(3) = (accel_intermediate(3)*I-m_C*g*norm(r_GC)*cos(th_GC_+th))/I;
        [accel, comparison(i)] = LBFcontrollerJerk(s0, param, accel_intermediate, deltaT, u, LfB_J, LgB_J, B_J, alpha, LfB, LgB, B);
        u = [m*accel(1); m*accel(2)+m*g; I*accel(3)+m_C*g*norm(r_GC)*cos(th_GC_+th)];
    elseif (controlMode == 8) % 8. LBF controller
        % states
        th = s0(5);
        r_GC = [s0(7); s0(8)]; % object frame

        % derived quantities
        m = m_G + m_C; % combined mass
        I = I_G + I_C + m_C*norm(r_GC)^2; % combined inertia
        th_GC_ = atan(r_GC(2)/r_GC(1)); % angle between r_GC and object frame x axis

        [~, ind] = min(abs(deltaT*i-tSave));
        sNominal = sSave(ind,:)';
        uNominal = uSave(i,:)';
        gains = [1; 1; 10]; % x y th
        alpha = 10e3;
        nu = 10e5;
        accel_impedence = impedenceController(s0,param,sNominal,uNominal,gains);
        accel_intermediate = accel_impedence-[0;g;0];
        accel_intermediate(3) = (accel_intermediate(3)*I-m_C*g*norm(r_GC)*cos(th_GC_+th))/I;
        [accel, comparison(i)] = LBFcontroller(s0, param, accel_intermediate, deltaT, u, LfB, LgB, B, alpha, nu);
        u = [m*accel(1); m*accel(2)+m*g; I*accel(3)+m_C*g*norm(r_GC)*cos(th_GC_+th)];
    elseif (controlMode == 9) % 9. Projection controller using 1 step optimization
        % states
        th = s0(5);
        r_GC = [s0(7); s0(8)]; % object frame

        % derived quantities
        m = m_G + m_C; % combined mass
        I = I_G + I_C + m_C*norm(r_GC)^2; % combined inertia
        th_GC_ = atan(r_GC(2)/r_GC(1)); % angle between r_GC and object frame x axis
        
        % impedence nominal controller fed into projection
        [~, ind] = min(abs(deltaT*i-tSave));
        sNominal = sSave(ind,:)';
        uNominal = uSave(i,:)';
        gains = [1; 1; 10]; % x y th
        accel_impedence = impedenceController(s0,param,sNominal,uNominal,gains);
        accel_intermediate = accel_impedence-[0;g;0];
        accel_intermediate(3) = (accel_intermediate(3)*I-m_C*g*norm(r_GC)*cos(th_GC_+th))/I;
        accel = controller(deltaT, 1, s0, param, fCone, vec, tol, jacob, accel_intermediate);
        u = [m*accel(1); m*accel(2)+m*g; I*accel(3)+m_C*g*norm(r_GC)*cos(th_GC_+th)];
    else % 10. MPC-DCBF controller
        % states
        th = s0(5);
        r_GC = [s0(7); s0(8)]; % object frame

        % derived quantities
        m = m_G + m_C; % combined mass
        I = I_G + I_C + m_C*norm(r_GC)^2; % combined inertia
        th_GC_ = atan(r_GC(2)/r_GC(1)); % angle between r_GC and object frame x axis

        % MPC parameters
        nHorizon = 10; % number of steps
        dt = deltaT*10; % time delta between each step

        %[~, ind] = min(abs(deltaT*i-tSave));
        %sNominal = sSave(ind,:)';
        %uNominal = uSave(i,:)';        
        gamma = ones(4,1).*0.8; % 0.4 works great for aggressive and 0.05 for conservative (or 0.1)
        %gains = [1; 1; 10]; % x y th
        %accel_impedence = impedenceController(s0,param,sNominal,uNominal,gains);
        %accel_intermediate = accel_impedence-[0;g;0];
        %accel_intermediate(3) = (accel_intermediate(3)*I-m_C*g*norm(r_GC)*cos(th_GC_+th))/I;
        if (i == 1)
            Bp_initial = [450; 450; 450; 450];
        end
        x_g = [5;0;0;0;0;0];
        %x_g = [2*sin(i*deltaT*3);0;2*(1-cos(i*deltaT*3));0;0;0]; % circle
        accel = controller_MPC_DCBF(dt, nHorizon, s0, param, Bp, Bp_initial, gamma, x_g);
        u = [m*accel(1); m*accel(2)+m*g; I*accel(3)+m_C*g*norm(r_GC)*cos(th_GC_+th)];
        Bp_initial = Bp([s0(1:6,1);accel]);
    end

    % Save control input
    if (saveController)
        uSave = [uSave; u'];
    end

    % simulation setup
    dynStick_ = @(t,s) dynStick(t,s,param,u);
    stickEvent_ = @(t,s) stickEvent(t,s,param,u,fCone,vec);
    odeOptionsStick = odeset('Events', stickEvent_);
    dynSlide_ = @(t,s) dynSlide(t,s,param,u);
    slideEvent_ = @(t,s) slideEvent(t,s,param,u,fCone,vec);
    odeOptionsSlide = odeset('Events', slideEvent_);

    % inner simulation time loop
    currentTime = (i-1)*deltaT;
    while(currentTime < i*deltaT)
        % checking initial constraint condition for dynamic grasping
        if (mode == 0)
            ie = find(checkSlip(worldToObject(0,s0,param,u),fCone,vec)>0,1);
            mode = modeSwitch(mode,ie);
            if (mode == -1)
                disp('Warning: initial action out of friction cone!')
                break;
            end
        end

        % run hybrid ode and add missing states back
        tOut = [];
        sOut = [];
        switch mode
            case 0
                s0 = s0(1:8);
                [tOut, sOut, ~, ~, ie] = ode45(dynStick_, [currentTime i*deltaT], s0, odeOptionsStick);
                sOut = [sOut zeros(size(sOut,1),4)];
            case 1
                s0 = s0(1:10);
                [tOut, sOut, ~, ~, ie] = ode45(dynSlide_, [currentTime i*deltaT], s0, odeOptionsSlide);
                sOut = [sOut zeros(size(sOut,1),2)];
            otherwise            
        end

        % record trajectory
        tTraj = [tTraj; tOut(1:end-1)];
        sTraj = [sTraj; sOut(1:end-1,:)];
        uTraj = [uTraj; repmat(u',length(tOut(1:end-1)),1)];
        mTraj = [mTraj; ones(size(sOut,1)-1,1).*mode];
            
        % reallocate initial conditions
        s0 = sOut(end,:)';
    
        % check if event function tripped
        mode = modeSwitch(mode,ie);
        currentTime = tOut(end);
        if (mode == -1)
            disp('Warning: state moved out of friction cone!')
            break;
        end
    end

    if (mode == 1 && ~trippedSlip)
        trippedSlip = true;
        dr_GC = [s0(9); s0(10)];
        [LfB1, LgB1, B1] = barrierFuncGenPlane(r_GC+dr_GC*deltaT*nHorizon, param, 1, fCone, vec, tol);
        [LfB2, LgB2, B2] = barrierFuncGenPlane(r_GC+dr_GC*deltaT*nHorizon, param, 2, fCone, vec, tol);
        [LfB3, LgB3, B3] = barrierFuncGenPlane(r_GC+dr_GC*deltaT*nHorizon, param, 3, fCone, vec, tol);
        [LfB4, LgB4, B4] = barrierFuncGenPlane(r_GC+dr_GC*deltaT*nHorizon, param, 4, fCone, vec, tol);
        Bp = @(s_u) [B1(s_u); B2(s_u); B3(s_u); B4(s_u)];
    end

    if (mode == -1)
        disp(['Error Time = ', num2str(currentTime)])
        break;
    else
        disp(['Time = ', num2str(currentTime)])
    end
end

if (saveController)
    tSave = tTraj;
    sSave = sTraj;
    save(fileName,'uSave','s0Save','tSave','sSave');
end

%% Animate Simulation
% save video options
videoFileName = [fileName(1:end-4) '-' date '-' num2str(controlMode) specialName '.avi'];
if (saveVideo)
    v = VideoWriter(videoFileName);
    v.Quality = 100;
end

% animation options
showAnimation = 1;
frictionAnimate = 0;
movingWindow = 0;

% plot options
plotbarrier = 1;

% interpolate dataset based on fps
fps = 120;
dt = 1/fps;
t = 0:dt:tTraj(end);
s = zeros(length(t),size(sTraj,2));
u = zeros(length(t),size(uTraj,2));
modes = -1.*ones(length(t),size(mTraj,2));
for i = 1:size(s,2)
    s(:,i) = interp1(tTraj, sTraj(:,i), t);
    %sS(:,i) = interp1(tSave, sSave(:,i), t);
end
for i = 1:size(u,2)
    for tt = 1:length(t)
        ind = find(tTraj >= t(tt),1);
        u(tt,i) = uTraj(ind,i);
    end
    %u(:,i) = interp1(tTraj, uTraj(:,i), t);
    %uS(:,i) = interp1(0:deltaT:10-deltaT, uSave(:,i), t);
end
for tt = 1:length(t)
    ind = find(tTraj >= t(tt),1);
    modes(tt) = mTraj(ind);
end
x_G = s(:,1);
dx_G = s(:,2);
y_G = s(:,3);
dy_G = s(:,4);
th = s(:,5);
dth = s(:,6);
r_GC = [s(:,7), s(:,8)];
CoR_radius = 0.005;

% initialize animation
if (showAnimation)
    figure(1000); clf;
    axis equal
end

if (saveVideo)
    open(v);
end
if (showAnimation)
    for i = 1:length(t)
        hold off
        hand = plot([-handW/2, handW/2, handW/2, -handW/2, -handW/2]+x_G(i), ...
            [handH/2, handH/2, -handH/2, -handH/2, handH/2]+y_G(i));
        hold on
        object = plot([-objectW/2, objectW/2, objectW/2, -objectW/2, -objectW/2]+x_G(i)+r_GC(i,1), ...
            [objectH/2, objectH/2, -objectH/2, -objectH/2, objectH/2]+y_G(i)+r_GC(i,2));
        CoR = rectangle('Position',[x_G(i)-CoR_radius/2 y_G(i)-CoR_radius/2 ...
            CoR_radius CoR_radius],'Curvature',[1 1]);
        rotate(hand,[0 0 1],rad2deg(th(i)),[x_G(i),y_G(i),0])
        rotate(object,[0 0 1],rad2deg(th(i)),[x_G(i),y_G(i),0])
        if (movingWindow)
            axis equal
            axis manual
            plot(x_G(1:i),y_G(1:i))
        else
            plot(x_G(1:i),y_G(1:i))
            xlim([min([x_G;y_G])-1 max([x_G;y_G])+1])
            ylim([min([x_G;y_G])-1 max([x_G;y_G])+1])
        end
        ylabel('y [m]')
        xlabel('x [m]')
        drawnow
        if (saveVideo)
            frame = getframe(gcf);
            writeVideo(v,frame);
        end
    end
end
if (saveVideo)
    close(v)
end

% Friction cone visualization
if (frictionAnimate)
    figure(1)
    for i = 1:length(t)
        hold off
        scatter3(fCone(1,:), fCone(2,:), fCone(3,:))
        hold on, scatter3(0,0,0)
        p = worldToObject(t(i),s(i,:)',param,u(i,:)');
        scatter3(p(1),p(2),p(3))
        quiver3([0 0 0 0]', [0 0 0 0]', [0 0 0 0]', fCone(1,:)', fCone(2,:)', fCone(3,:)',1)
        pause(0.01)
    end
else
    figure(1)
    %scatter3(fCone(1,:), fCone(2,:), fCone(3,:))
    hold on%, scatter3(0,0,0)
    %quiver3([0 0 0 0]', [0 0 0 0]', [0 0 0 0]', fCone(1,:)', fCone(2,:)', fCone(3,:)',1)
    patch([0; fCone(1,1); fCone(1,2)],[0;fCone(2,1);fCone(2,2)],[0;fCone(3,1);fCone(3,2)],[0.5,0,0.5],'FaceAlpha',0.1)
    patch([0; fCone(1,2); fCone(1,3)],[0;fCone(2,2);fCone(2,3)],[0;fCone(3,2);fCone(3,3)],[0.5,0,0.5],'FaceAlpha',0.1)
    patch([0; fCone(1,3); fCone(1,4)],[0;fCone(2,3);fCone(2,4)],[0;fCone(3,3);fCone(3,4)],[0.5,0,0.5],'FaceAlpha',0.1)
    patch([0; fCone(1,4); fCone(1,1)],[0;fCone(2,4);fCone(2,1)],[0;fCone(3,4);fCone(3,1)],[0.5,0,0.5],'FaceAlpha',0.1)
    p = zeros(3,length(t));
    for i = 1:length(t)
        p(:,i) = worldToObject(t(i),s(i,:)',param,u(i,:)');
    end
    %for i = 1:length(t)
    %    pSave(:,i) = worldToObject(t(i),sS(i,:)',param,uS(i,:));
    %end
    scatter3(p(1,1),p(2,1),p(3,1),'k');
    plot3(p(1,:),p(2,:),p(3,:),'r')
    %hold on
    %scatter3(pSave(1,1),pSave(2,1),pSave(3,1),'k')
    %plot3(pSave(1,:),pSave(2,:),pSave(3,:),'b')
    %legend({'','','','','','\gamma = 0.1','','\gamma = 0.4'})
    xlabel('$\ddot{x}_{c} [m/s^2]$','interpreter','latex')
    ylabel('$\ddot{y}_{c} [m/s^2]$','interpreter','latex')
    zlabel('$\ddot{\theta}_{c} [rad/s^2]$','interpreter','latex')
    title('Gravito-Inertial Wrench Constraints')
    if (plotbarrier)
        mCH = maxNormalForce/m_C;
        fConeImplicit = @(p1,p2,p3) - p1.^2/a - p3.^2/b + p2.^2;
        fPlaneImplicit = @(p1,p2,p3,i) vec(1,i)*(p1-fCone(1,i)) + vec(2,i)*(p2-fCone(2,i)) + ...
            vec(3,i)*(p3-fCone(3,i));
        if (controlMode == 7 || controlMode == 8)
            fimplicit3(fConeImplicit,[-mCH, mCH, 0, mCH, -mCH, mCH],'EdgeColor','None','FaceAlpha',0.1);
        end
        %fimplicit3(@(p1,p2,p3) fPlaneImplicit(p1,p2,p3,1), [-mCH, mCH, 0, mCH, -mCH, mCH],'EdgeColor','None','FaceAlpha',0.1);
        %fimplicit3(@(p1,p2,p3) fPlaneImplicit(p1,p2,p3,2), [-mCH, mCH, 0, mCH, -mCH, mCH],'EdgeColor','None','FaceAlpha',0.1);
        %fimplicit3(@(p1,p2,p3) fPlaneImplicit(p1,p2,p3,3), [-mCH, mCH, 0, mCH, -mCH, mCH],'EdgeColor','None','FaceAlpha',0.1);
        %fimplicit3(@(p1,p2,p3) fPlaneImplicit(p1,p2,p3,4), [-mCH, mCH, 0, mCH, -mCH, mCH],'EdgeColor','None','FaceAlpha',0.1);
        % CBF
        if (controlMode == 7 || controlMode == 8)
            figure(100)
            Bplot = fConeImplicit(p(1,:),p(2,:),p(3,:));
            plot(t,Bplot)
            title('CBF')
            xlabel('time [s]')
            ylabel('B(x,u)')
            hold on, yyaxis right
            tT = linspace(0,deltaT*length(comparison),length(comparison));
            ylabel('B approximation')
            Bplotapprox = cumsum(comparison);
            plot(tT, Bplotapprox)
        elseif (controlMode == 10)
            figure(100)
            Bplot1 = -fPlaneImplicit(p(1,:),p(2,:),p(3,:),1);
            Bplot2 = -fPlaneImplicit(p(1,:),p(2,:),p(3,:),2);
            Bplot3 = -fPlaneImplicit(p(1,:),p(2,:),p(3,:),3);
            Bplot4 = -fPlaneImplicit(p(1,:),p(2,:),p(3,:),4);
            plot(t,Bplot1), hold on
            plot(t,Bplot2)
            plot(t,Bplot3)
            plot(t,Bplot4)
            title('CBF')
            xlabel('time [s]')
            ylabel('h(s,u)')
        end
    end
end

%% Plots
% Trajectory
figure(2)
if (any(stabilizingControl == controlMode))
    plot(x_G,y_G)
    hold on,
    plot(sSave(:,1),sSave(:,3))
    xlabel('x [m]')
    ylabel('y [m]')
    legend({'online execution', 'nominal trajectory'})
    title('Spatial Trajectory')
else
    plot(x_G,y_G)
    xlabel('x [m]')
    ylabel('y [m]')
    title('Spatial Trajectory')
end

% Orientation
figure(3)
if (any(stabilizingControl == controlMode))
    plot(t,th)
    hold on,
    plot(tSave,sSave(:,5))
    xlabel('time [s]')
    ylabel('angle [rad]')
    legend({'online execution', 'nominal trajectory'})
    title('Orientation')
else
    plot(t,th)
    xlabel('time [s]')
    ylabel('angle [rad]')
    title('Orientation')
end


% Control input
figure(4)
plot(t,u)
xlabel('time [s]')
ylabel('control effort')
legend({'u_x','u_y','u_{th}'})
title('Control Inputs')