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
wE = 0;

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
        if (thetaCheck-pi/2) <= phi(i,t) && phi(i,t) <= (thetaCheck+pi/2)
            rho(i,t) = sqrt((X(t)-Xs(i,t))^2 + (Y(t)-Ys(i,t))^2);
            rhoDot(i,t) = ((X(t)-Xs(i,t))*(XD(t)-XDs(i,t)) + (Y(t)-Ys(i,t))*(YD(t)-YDs(i,t)))...
                / rho(i,t);
        else
            rho(i,t) = nan;
            rhoDot(i,t) = nan;
        end
    end
end
%test plots for measurement
fig=figure; hold on; grid on; grid minor
for i=1:12
    plot(T,rho(i,:)) 
end
legend('S1','S2','S3','S4','S5','S6','S7','S8','S9','S10','S11','S12')
title('\rho measurements by station'); xlabel('Time [s]'); ylabel('\rho [km]')
saveas(fig,'ASEN5044_HW8_P2_rhoModel.png','png');

fig=figure; hold on; grid on; grid minor
for i=1:12
    plot(T,rhoDot(i,:))
end
legend('S1','S2','S3','S4','S5','S6','S7','S8','S9','S10','S11','S12')
title('\rhoDot measurements by station'); xlabel('Time [s]'); ylabel('\rhoDot [km/s]')
saveas(fig,'ASEN5044_HW8_P2_rhoDotModel.png','png');

fig=figure; hold on; grid on; grid minor
for i=1:12
    plot(T,phi(i,:))
end
legend('S1','S2','S3','S4','S5','S6','S7','S8','S9','S10','S11','S12')
title('\phi measurements by station'); xlabel('Time [s]'); ylabel('\phi [rad]')
saveas(fig,'ASEN5044_HW8_P2_phi.png','png');

dt = 10;

% Nominalize the A matrix
delta_x = [0 0 0 0]';

dx = [0.01, 0.01, 0.01, 0.01]';
dx_lin = [];
for i = 1:length(t_vec)
    dx_lin = horzcat(dx_lin, dx(:,i));
    %based on current position, calculate F
    F = F_variant(dx_lin(1,i),dx_lin(3,i));
    
    %using that, figure out next dx
    dx(:,i+1) = F*dx(:,i);
end
dx_lin = dx_lin';

%====== c ======%
% One Full Orbit Period
% Simulate the Linearized DT dynamics


fprintf('Plotting Satellite for 1 Orbit:\n');
fprintf("x0 = [6678, 0, 0, 7.7258]'");

fig = figure('visible','on');
set(fig,'Position',[100 100 900 600]);
sgtitle('Satellite State');

subplot(2,2,1); hold on; grid on; grid minor;
plot(T,S(:,1),'b-','LineWidth',1.5);
plot(t_vec,S(:,1)+dx_lin(:,1),'r--','LineWidth',1.5);
ylim([-7000 7000]);
xlabel('time [sec]');
ylabel('X position [km]');
legend('ODE45 DT Nonlinear Simulation (dt=10s)','DT Linearized Simulation (dt=10s)');
xlim([0 P]);

subplot(2,2,2); hold on; grid on; grid minor;
plot(T,S(:,3),'b-','LineWidth',1.5);
plot(t_vec,S(:,3)+dx_lin(:,3),'r--','LineWidth',1.5);
ylim([-7000 7000]);
xlabel('time [sec]');
ylabel('Y position [km]');
xlim([0 P]);

subplot(2,2,3); hold on; grid on; grid minor;
plot(T,S(:,2),'b-','LineWidth',1.5);
plot(t_vec,S(:,2)+dx_lin(:,2),'r--','LineWidth',1.5);
ylim([-10 10]);
xlabel('time [sec]');
ylabel('X velocity [km/s]');
xlim([0 P]);

subplot(2,2,4); hold on; grid on; grid minor;
plot(T,S(:,4),'b-','LineWidth',1.5);
plot(t_vec,S(:,4)+dx_lin(:,4),'r--','LineWidth',1.5);
ylim([-10 10]);
xlabel('time [sec]');
ylabel('Y velocity [km/s]');
xlim([0 P]);

saveas(fig,'ASEN5044_HW8_P2_sim.png','png');


fig = figure('visible','on');
set(fig,'Position',[100 100 900 600]);
sgtitle('Satellite State - Linearized');

subplot(2,2,1); hold on; grid on; grid minor;
plot(t_vec,dx_lin(:,1),'r--','LineWidth',1.5);
xlabel('time [sec]');
ylabel('dX position [km]');

subplot(2,2,2); hold on; grid on; grid minor;
plot(t_vec,dx_lin(:,3),'r--','LineWidth',1.5);
xlabel('time [sec]');
ylabel('dY position [km]');
xlim([0 P]);

subplot(2,2,3); hold on; grid on; grid minor;
plot(t_vec,dx_lin(:,2),'r--','LineWidth',1.5);
xlabel('time [sec]');
ylabel('dX velocity [km/s]');
xlim([0 P]);

subplot(2,2,4); hold on; grid on; grid minor;
plot(t_vec,dx_lin(:,4),'r--','LineWidth',1.5);
xlabel('time [sec]');
ylabel('dY velocity [km/s]');
xlim([0 P]);

saveas(fig,'ASEN5044_HW8_P2_lin.png','png');



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

H = [2*sqrt(X-Xs), 0, 2*sqrt(Y-Ys), 0;
        ((Xdot-Xsdot)/(sqrt((X-Xs)^2+(Y-Ys)^2)))-((2*(X-Xs)*((X-Xs)*(Xdot-Xsdot)+(Y-Ys)*(Ydot-Ysdot)))/(sqrt((X-Xs)^2+(Y-Ys)^2))),...
        ((X-Xs)/(sqrt((X-Xs)^2+(X-Xs)^2))),...
        -((2*(Y-Ys)*((X-Xs)*(Xdot-Xsdot)+(Y-Ys)*(Ydot-Ysdot)))/(sqrt((X-Xs)^2+(Y-Ys)^2)))+((Ydot-Ysdot)/(sqrt((X-Xs)^2+(Y-Ys)^2))),...
        ((Y-Ys)/(sqrt((X-Xs)^2+(Y-Ys)^2)));...
        ((Y-Ys)/((X-Xs)^2+(Y-Ys)^2)),0,...
        ((-X+Xs)/((X-Xs)^2+(Y-Ys)^2)),0];
        
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
