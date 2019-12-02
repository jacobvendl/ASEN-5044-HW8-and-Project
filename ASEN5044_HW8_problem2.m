%==========================================================================
% Jake Vendl | Jack Toland
% ASEN 5044
% Homework 8
% 12/3/2019
%==========================================================================
close all; clear all; clc

mu = 398600;        % km^3/s^2
r0 = 6678;          % km
rE = 6378;          % km
wE = 2*pi/86400;    % rad/s

x0 = [6678, 0, 0, r0*sqrt(mu/r0^3)]';
P = 2*pi*sqrt(r0^3/mu);

t_vec = 0:10:P;

%====== a ======%

r0_nom = 6678;
x1_nom = 6678;
x3_nom = 0;
Atil = [0, 1, 0, 0;
        (-mu/r0_nom^3)+((3*mu*x1_nom^2)/r0_nom^5), 0, ((3*mu*x1_nom*x3_nom)/r0_nom^5), 0;
        0, 0, 0, 1;
        ((3*mu*x1_nom*x3_nom)/r0_nom^5), 0, (-mu/r0^3)+((3*mu*x3_nom^3)/r0_nom^5), 0];
    
Btil = [0, 0;
        1, 0;
        0, 0;
        0, 1];
    
%Ctil = [2*sqrt(r0-xs), 0, 2*sqrt(-ys), 0;
        


%====== b ======%
s0 = x0;
opts = odeset('RelTol',1e-12,'AbsTol',1e-12);
[T, S] = ode45(@(t,s)orbit_prop_func(t,s),t_vec,s0,opts);

X=S(:,1); Y=S(:,3); XD=S(:,2); YD=S(:,4);
Xs = zeros(12,length(T));
Ys = zeros(12,length(T));
XDs = zeros(12,length(T));
YDs = zeros(12,length(T));
rho = zeros(12,length(T));
rhoDot = zeros(12,length(T));
phi = zeros(12,length(T));
%now simulate the measurements for all time
for i=1:12 %stations
    theta = (i-1)*pi/6;
    for t=1:length(T) %loop through one orbit period
        currentTime = T(t);
        
        %find station position and velocity
        Xs(i,t) = rE*cos(wE*currentTime + theta);
        Ys(i,t) = rE*sin(wE*currentTime + theta);
        XDs(i,t) = -rE*wE*sin(wE*currentTime + theta);
        YDs(i,t) = rE*wE*cos(wE*currentTime + theta);
        
        %peform check at given time to see if s/c is visible
        phi(i,t) = atan2((Y(t)-Ys(i,t)),(X(t)-Xs(i,t)));
        thetaCheck = atan2(Ys(i,t),Xs(i,t));
        if (thetaCheck-pi/2) > (thetaCheck+pi/2)
            upperBound = thetaCheck-pi/2;
            lowerBound = thetaCheck+pi/2;
        else
            upperBound = thetaCheck+pi/2;
            lowerBound = thetaCheck-pi/2;
        end
        if (lowerBound <= phi(i,t) && phi(i,t) <= upperBound) ...
            || (lowerBound-2*pi <= phi(i,t) && phi(i,t)<=upperBound-2*pi)...
            || (lowerBound+2*pi <= phi(i,t) && phi(i,t)<=upperBound+2*pi)
            
            rho(i,t) = sqrt((X(t)-Xs(i,t))^2 + (Y(t)-Ys(i,t))^2);
            rhoDot(i,t) = ((X(t)-Xs(i,t))*(XD(t)-XDs(i,t)) + (Y(t)-Ys(i,t))*(YD(t)-YDs(i,t)))...
                / rho(i,t);
        else
            rho(i,t) = nan;
            rhoDot(i,t) = nan;
        end
    end
end

dt = 10;

% Nominalize the A matrix. starting with initial perturbance
dx = [0.1, 0.01, 0.1, 0.01]';
dx0 = dx;
dx_lin = [];
dy_lin = zeros(36,length(T));
for t = 1:length(t_vec)
    dx_lin = horzcat(dx_lin, dx(:,t));
    %based on current position, calculate F
    F = F_variant(dx_lin(1,t),dx_lin(3,t));
    
    %using that, figure out next dx
    dx(:,t+1) = F*dx(:,t);
    
    %now find linearized measurement
    for i=1:12
        if ~isnan(rho(i,t))
            H = H_variant(X(t),XD(t),Y(t),YD(t),Xs(i,t),XDs(i,t),Ys(i,t),YDs(i,t));
            dy(:,t) = H*dx(:,t);
            dy_lin(3*i-2:3*i,t) = dy(:,t);
        else
            dy_lin(3*i-2:3*i,t) = [nan nan nan]';
        end
    end
    
end
dx_lin = dx_lin';

%====== c ======%
% One Full Orbit Period
% Simulate the Linearized DT dynamics

s0 = x0 + dx0;
opts = odeset('RelTol',1e-12,'AbsTol',1e-12);
[T_per, S_per] = ode45(@(t,s)orbit_prop_func(t,s),t_vec,s0,opts);


fprintf('Plotting Satellite for 1 Orbit:\n');
fprintf("x0 = [6678, 0, 0, 7.7258]'");

fig = figure('visible','on');
set(fig,'Position',[100 100 900 600]);
sgtitle(sprintf('Satellite State Model - Nonlinear and Linear Simulation \n dx0 = [%.2fkm %.2fkm/s %.2fkm %.2fkm/s]',dx(1,1),dx(2,1),dx(3,1),dx(4,1)));

subplot(2,2,1); hold on; grid on; grid minor;
plot(T_per,S_per(:,1),'b-','LineWidth',1.2);
plot(t_vec,S(:,1)+dx_lin(:,1),'r--','LineWidth',1.2);
ylim([-7000 7000]);
xlabel('time [sec]');
ylabel('X position [km]');
legend('ODE45 DT Nonlinear Simulation (dt=10s)','DT Linearized Simulation (dt=10s)');
xlim([0 P]);

subplot(2,2,2); hold on; grid on; grid minor;
plot(T_per,S_per(:,3),'b-','LineWidth',1.5);
plot(t_vec,S(:,3)+dx_lin(:,3),'r--','LineWidth',1.5);
ylim([-7000 7000]);
xlabel('time [sec]');
ylabel('Y position [km]');
xlim([0 P]);

subplot(2,2,3); hold on; grid on; grid minor;
plot(T_per,S_per(:,2),'b-','LineWidth',1.5);
plot(t_vec,S(:,2)+dx_lin(:,2),'r--','LineWidth',1.5);
ylim([-10 10]);
xlabel('time [sec]');
ylabel('X velocity [km/s]');
xlim([0 P]);

subplot(2,2,4); hold on; grid on; grid minor;
plot(T_per,S_per(:,4),'b-','LineWidth',1.5);
plot(t_vec,S(:,4)+dx_lin(:,4),'r--','LineWidth',1.5);
ylim([-10 10]);
xlabel('time [sec]');
ylabel('Y velocity [km/s]');
xlim([0 P]);

saveas(fig,'ASEN5044_HW8_P2_sim.png','png');


fig = figure('visible','on');
set(fig,'Position',[100 100 900 600]);
sgtitle(sprintf('Satellite State Model - Linearized Deviations \n dx0 = [%.2fkm %.2fkm/s %.2fkm %.2fkm/s]',dx(1,1),dx(2,1),dx(3,1),dx(4,1)));

subplot(2,2,1); hold on; grid on; grid minor;
plot(t_vec,S_per(:,1)-S(:,1),'b-','LineWidth',1.5);
plot(t_vec,dx_lin(:,1),'r--','LineWidth',1.5);
xlabel('time [sec]');
ylabel('\deltaX position [km]');
legend('ODE45 DT Nonlinear Simulation (dt=10s)','DT Linearized Simulation (dt=10s)','Location','Southwest');


subplot(2,2,2); hold on; grid on; grid minor;
plot(t_vec,S_per(:,3)-S(:,3),'b-','LineWidth',1.5);
plot(t_vec,dx_lin(:,3),'r--','LineWidth',1.5);
xlabel('time [sec]');
ylabel('\deltaY position [km]');
xlim([0 P]);

subplot(2,2,3); hold on; grid on; grid minor;
plot(t_vec,S_per(:,2)-S(:,2),'b-','LineWidth',1.5);
plot(t_vec,dx_lin(:,2),'r--','LineWidth',1.5);
xlabel('time [sec]');
ylabel('\deltaVX [km/s]');
xlim([0 P]);

subplot(2,2,4); hold on; grid on; grid minor;
plot(t_vec,S_per(:,4)-S(:,4),'b-','LineWidth',1.5);
plot(t_vec,dx_lin(:,4),'r--','LineWidth',1.5);
xlabel('time [sec]');
ylabel('\deltaVY [km/s]');
xlim([0 P]);

saveas(fig,'ASEN5044_HW8_P2_lin.png','png');


fig = figure('visible','on'); 
set(fig,'Position',[100 100 900 600]);
sgtitle(sprintf('Satellite Measurement Model - Nonlinear and Linear Simulation \n dx0 = [%.2f %.2f %.2f %.2f]',dx(1,1),dx(2,1),dx(3,1),dx(4,1)));
subplot(3,1,2); hold on; grid on; grid minor;
for i=1:12
    plot(T,rhoDot(i,:),'-','LineWidth',1.5)
    plot(T,rhoDot(i,:) + dy_lin(3*i-1,:),'k--','LineWidth',1.5)
end
xlabel('Time [s]'); ylabel('\rhoDot^i [km/s]')

subplot(3,1,3); hold on; grid on; grid minor;
for i=1:12
    plot(T,phi(i,:),'-','LineWidth',1.5)
    plot(T,phi(i,:) + dy_lin(3*i,:),'k--','LineWidth',1.5)
end
subplot(3,1,1); hold on; grid on; grid minor;
for i=1:12
    plot(T,rho(i,:),'-','LineWidth',1.5)
    plot(T,rho(i,:) + dy_lin(3*i-2,:),'k--','LineWidth',1.5)
end
xlabel('Time [s]'); ylabel('\rho^i [km]')
legend('S1','S1_{lin}','S2','S2_{lin}','S3','S3_{lin}','S4','S4_{lin}','S5','S5_{lin}'...
    ,'S6','S6_{lin}','S7','S7_{lin}','S8','S8_{lin}','S9','S9_{lin}'...
    ,'S10','S10_{lin}','S11','S11_{lin}','S12','S12_{lin}')
xlabel('Time [s]'); ylabel('\phi^i [rad]')

saveas(fig,'ASEN5044_HW8_P2_MeasurementModel.png','png');


fig = figure('visible','on'); 
set(fig,'Position',[100 100 900 600]);
sgtitle(sprintf('Satellite Measurement Model - Linearized Deviations \n dx0 = [%.2f %.2f %.2f %.2f]',dx(1,1),dx(2,1),dx(3,1),dx(4,1)));
subplot(3,1,2); hold on; grid on; grid minor;
for i=1:12
    plot(T,dy_lin(3*i-1,:),'--','LineWidth',1.5)
end
xlabel('Time [s]'); ylabel('\delta\rhoDot^i [km/s]')
subplot(3,1,3); hold on; grid on; grid minor;
for i=1:12
    plot(T,dy_lin(3*i,:),'--','LineWidth',1.5)
end
subplot(3,1,1); hold on; grid on; grid minor;
for i=1:12
    plot(T,dy_lin(3*i-2,:),'--','LineWidth',1.5)
end
xlabel('Time [s]'); ylabel('\delta\rho^i [km]')
legend('S1_{lin}','S2_{lin}','S3_{lin}','S4_{lin}','S5_{lin}'...
    ,'S6_{lin}','S7_{lin}','S8_{lin}','S9_{lin}'...
    ,'S10_{lin}','S11_{lin}','S12_{lin}')
xlabel('Time [s]'); ylabel('\delta\phi^i [rad]')

saveas(fig,'ASEN5044_HW8_P2_LinearizedMeasurementModel.png','png');


%make movie simulation
% fig=figure; hold on; grid on; grid minor;
% for t=1:length(T)
%     fig; hold on; grid on; grid minor;
%     plot(X(t),Y(t),'b*')
%     for i=1:12
%         plot(Xs(i,t),Ys(i,t),'k.')
%         if ~isnan(rho(i,t))
%             line([Xs(i,t) X(t)],[Ys(i,t) Y(t)])
%         end
%     end
%     
%     axis equal; axis([-8000 8000 -8000 8000]); 
%     xlabel('x');ylabel('y');title(sprintf('simulation, t=%.0fs',T(t)))
%     M(t) = getframe(fig);
%     clf(fig);
% end
% video = VideoWriter('visualization', 'Uncompressed AVI');
% video.FrameRate = 60;
% open(video)
% writeVideo(video, M);
% close(video);

function [ F ] = F_variant(X,Y)

mu = 398600;        % km^3/s^2
r0_nom = 6678;          % km
dt = 10;

F = expm(dt*[0, 1, 0, 0;
        (-mu*(r0_nom)^(-3))+(3*mu*X^2*r0_nom^(-5)), 0, 3*mu*X*Y*r0_nom^(-5), 0;
        0, 0, 0, 1;
        (3*mu*X*Y)*r0_nom^(-5), 0, (-mu*r0_nom^(-3))+(3*mu*Y^2*r0_nom^(-5)), 0]);
end


function [ H ] = H_variant(X,Xdot,Y,Ydot,Xs,Xsdot,Ys,Ysdot)
%initialize H
H = zeros(3,4);
    
%first row
H(1,1) = (X-Xs)/sqrt((X-Xs)^2+(Y-Ys)^2);
H(1,2) = 0;
H(1,3) = (Y-Ys)/sqrt((X-Xs)^2+(Y-Ys)^2);
H(1,4) = 0;

%second row
H(2,1) = (Xdot-Xsdot)/sqrt((X-Xs)^2+(Y-Ys)^2) - (X-Xs)*((X-Xs)*(Xdot-Xsdot)+(Y-Ys)*(Ydot-Ysdot)) / ((X-Xs)^2+(Y-Ys)^2)^(3/2);
H(2,3) = (Ydot-Ysdot)/sqrt((X-Xs)^2+(Y-Ys)^2) - (Y-Ys)*((X-Xs)*(Xdot-Xsdot)+(Y-Ys)*(Ydot-Ysdot)) / ((X-Xs)^2+(Y-Ys)^2)^(3/2);
H(2,2) = (X-Xs)/sqrt((X-Xs)^2+(Y-Ys)^2);
H(2,4) = (Y-Ys)/sqrt((X-Xs)^2+(Y-Ys)^2);
    
%third row
H(3,1) = ((Y-Ys)/((X-Xs)^2+(Y-Ys)^2));
H(3,2) = 0;
H(3,3) = ((-X+Xs)/((X-Xs)^2+(Y-Ys)^2));
H(3,4) = 0;
end

function [ ds ] = orbit_prop_func(t,s)

mu = 398600; 

x = s(1);
y = s(3);
r = sqrt(x^2+y^2);

xdot = s(2);
ydot = s(4);

xddot = -mu/r^3 * x;
yddot = -mu/r^3 * y;

ds = [xdot, xddot, ydot, yddot]';
end






