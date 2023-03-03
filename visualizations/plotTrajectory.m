function [] = plotTrajectory(s0, ss0, xp, yp, psolve, fCone, vec, r_GC, param, showSnapshots, movingWindow, ...
    objectH, objectW, handH, handW, xForceStretch, showForces, mu)

%% Spline
dxp = fnder(xp,1);
ddxp = fnder(xp,2);
dyp = fnder(yp,1);
ddyp = fnder(yp,2);

%% TOPP Results
P = reshape(psolve,numel(psolve)/5,5);
ss = ss0;
dss = P(:,1);
ddss = P(:,2);
ths = P(:,3);
dths = P(:,4);
ddths = P(:,5);

figure, plot(P(:,3))
title ('Orientation')
ylabel('angle [rad]')
figure, plot(P(:,1)), hold on, plot(P(:,2))
legend({'$\dot{s}$','$\ddot{s}$'},'Interpreter','latex')
title('$\dot{s}$ and $\ddot{s}$','Interpreter','latex')

%% Friction Cone

figure(100)
hold on
patch([0; fCone(1,1); fCone(1,2)],[0;fCone(2,1);fCone(2,2)],[0;fCone(3,1);fCone(3,2)],[0.5,0,0.5],'FaceAlpha',0.1)
patch([0; fCone(1,2); fCone(1,3)],[0;fCone(2,2);fCone(2,3)],[0;fCone(3,2);fCone(3,3)],[0.5,0,0.5],'FaceAlpha',0.1)
patch([0; fCone(1,3); fCone(1,4)],[0;fCone(2,3);fCone(2,4)],[0;fCone(3,3);fCone(3,4)],[0.5,0,0.5],'FaceAlpha',0.1)
patch([0; fCone(1,4); fCone(1,1)],[0;fCone(2,4);fCone(2,1)],[0;fCone(3,4);fCone(3,1)],[0.5,0,0.5],'FaceAlpha',0.1)
p = zeros(3,length(ss));
for i = 1:length(ss)
    ds = P(i,1);
    dds = P(i,2);
    th = P(i,3);
    dth = P(i,4);
    ddth = P(i,5);
    ax_G = fnval(ddxp,ss(i))*ds^2 + fnval(dxp,ss(i))*dds;
    ay_G = fnval(ddyp,ss(i))*ds^2 + fnval(dyp,ss(i))*dds;
    g = param(9);
    R_GC = norm(r_GC);
    th_GC = angle(r_GC(1) + 1i*r_GC(2));
    p(:,i) = [cos(th) sin(th) -R_GC*sin(th_GC); ...
              -sin(th) cos(th) R_GC*cos(th_GC); ...
              0 0 1]*[ax_G; ay_G; ddth] + ...
              [g*sin(th)-R_GC*dth^2*cos(th_GC); ...
              g*cos(th)-R_GC*dth^2*sin(th_GC); 0];
end
%for i = 1:length(t)
%    pSave(:,i) = worldToObject(t(i),sS(i,:)',param,uS(i,:));
%end
scatter3(p(1,1),p(2,1),p(3,1),'k');
plot3(p(1,:),p(2,:),p(3,:),'b')
%hold on
%scatter3(pSave(1,1),pSave(2,1),pSave(3,1),'k')
%plot3(pSave(1,:),pSave(2,:),pSave(3,:),'b')
%legend({'','','','','','\gamma = 0.1','','\gamma = 0.4'})
xlabel('$\ddot{x}_{c} [m/s^2]$','interpreter','latex')
ylabel('$\ddot{y}_{c} [m/s^2]$','interpreter','latex')
zlabel('$\ddot{\theta}_{c} [rad/s^2]$','interpreter','latex')
title('Gravito-Inertial Wrench Constraints')

% forward time scaling
dT = 2.*(ss(2:end)-ss(1:end-1))./(dss(2:end)+dss(1:end-1));
tss = cumsum([0; dT]);
dt = min(dT)/10;
dx = fnval(dxp,ss).*dss;
dy = fnval(dyp,ss).*dss;
dth = dths;
t = 0:dt:tss(end);
dx = interp1(tss,dx,t);
dy = interp1(tss,dy,t);
dth = interp1(tss,dth,t);
x_G = cumsum(dx.*dt)' + s0(1);
y_G = cumsum(dy.*dt)' + s0(3);
th = cumsum(dth.*dt)' + s0(5);
forces = individualForces(vec,p);
forces = forces./max(forces).*0.9*objectH;
fLx = forces(1,:); fLy = forces(2,:);
fRx = forces(3,:); fRy = forces(4,:);
fLx = interp1(tss,fLx,t); fLy = interp1(tss,fLy,t);
fRx = interp1(tss,fRx,t); fRy = interp1(tss,fRy,t);
fLx = fLx*xForceStretch; fRx = fRx*xForceStretch;

% initialize animation
CoR_radius = 0.005;
figure(1002); clf;
axis equal

% show animation
for i = 1:1+29*showSnapshots:length(t)
    if (showSnapshots)
        hold on
    else
        hold off
    end
    hand = plot([-handW/2, handW/2, handW/2, -handW/2, -handW/2]+x_G(i), ...
        [handH/2, handH/2, -handH/2, -handH/2, handH/2]+y_G(i), 'color', [0    0.4470    0.7410]);
    hold on
    object = plot([-objectW/2, objectW/2, objectW/2, -objectW/2, -objectW/2]+x_G(i)+r_GC(1), ...
        [objectH/2, objectH/2, -objectH/2, -objectH/2, objectH/2]+y_G(i)+r_GC(2), 'color', 'black'); % [0.8500    0.3250    0.0980]
    CoR = rectangle('Position',[x_G(i)-CoR_radius/2 y_G(i)-CoR_radius/2 ...
        CoR_radius CoR_radius],'Curvature',[1 1]);
    rotate(hand,[0 0 1],rad2deg(th(i)),[x_G(i),y_G(i),0])
    rotate(object,[0 0 1],rad2deg(th(i)),[x_G(i),y_G(i),0])
    if (showForces)
        forceL = plot([-objectW/2 -objectW/2+fLx(i)]+x_G(i)+r_GC(1), ...
            [-objectH/2 -objectH/2+fLy(i)]+y_G(i)+r_GC(2), 'color', 'red');
        forceR = plot([objectW/2 objectW/2+fRx(i)]+x_G(i)+r_GC(1), ...
            [-objectH/2 -objectH/2+fRy(i)]+y_G(i)+r_GC(2), 'color', 'red');
        frictionConeL = plot([-objectW/2 -objectW/2-objectH/2*mu*xForceStretch -objectW/2 -objectW/2+objectH/2*mu*xForceStretch]+x_G(i)+r_GC(1), ...
            [-objectH/2 0 -objectH/2 0]+y_G(i)+r_GC(2), '--', 'color', 'green');
        frictionConeR = plot([objectW/2 objectW/2-objectH/2*mu*xForceStretch objectW/2 objectW/2+objectH/2*mu*xForceStretch]+x_G(i)+r_GC(1), ...
            [-objectH/2 0 -objectH/2 0]+y_G(i)+r_GC(2), '--', 'color', 'green');
        rotate(forceL,[0 0 1],rad2deg(th(i)),[x_G(i),y_G(i),0])
        rotate(forceR,[0 0 1],rad2deg(th(i)),[x_G(i),y_G(i),0])
        rotate(frictionConeL,[0 0 1],rad2deg(th(i)),[x_G(i),y_G(i),0])
        rotate(frictionConeR,[0 0 1],rad2deg(th(i)),[x_G(i),y_G(i),0])
    end
    if (movingWindow)
        axis equal
        axis manual
        plot(x_G(1:i),y_G(1:i))
    else
        plot(x_G(1:i),y_G(1:i),'--')
        xlim([min([x_G;y_G])-1 max([x_G;y_G])+1])
        ylim([min([x_G;y_G])-1 max([x_G;y_G])+1])
    end
    ylabel('y [m]')
    xlabel('x [m]')
    drawnow
end


end

