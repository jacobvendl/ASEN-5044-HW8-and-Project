%==========================================================================
% Jake Vendl | Jack Toland
% ASEN 5044
% Homework 8
% 12/3/2019
%==========================================================================
close all; clear all; clc

mu = 398600;        % km^3/s^2
r0 = 6678;          % km

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

dt = 10;

% Nominalize the A matrix


updatex = x0;

x_sim = [];
for i = 1:length(t_vec)
        x_sim = horzcat(x_sim, updatex);
        F = F_variant(S(i,1),S(i,3));
        updatex = F*updatex;
end
x_sim = x_sim';



%====== c ======%
% One Full Orbit Period
% Simulate the Linearized DT dynamics





fprintf('Plotting Satellite for 1 Orbit:\n');
fprintf("x0 = [6678, 0, 0, 7.7258]' + [1, 1, 0.1, 0.1]'");

fig = figure('visible','on');
set(fig,'Position',[100 100 900 600]);
sgtitle('Satellite State');

subplot(2,2,1); hold on; grid on; grid minor;
plot(T,S(:,1),'b-','LineWidth',1.5);
plot(t_vec,x_sim(:,1),'r--','LineWidth',1.5);
ylim([-7000 7000]);
xlabel('time [sec]');
ylabel('X position [km]');
xlim([0 P]);

subplot(2,2,2); hold on; grid on; grid minor;
plot(T,S(:,3),'b-','LineWidth',1.5);
plot(t_vec,x_sim(:,3),'r--','LineWidth',1.5);
ylim([-7000 7000]);
xlabel('time [sec]');
ylabel('Y position [km]');
xlim([0 P]);

subplot(2,2,3); hold on; grid on; grid minor;
plot(T,S(:,2),'b-','LineWidth',1.5);
plot(t_vec,x_sim(:,2),'r--','LineWidth',1.5);
ylim([-10 10]);
xlabel('time [sec]');
ylabel('X velocity [km/s]');
xlim([0 P]);

subplot(2,2,4); hold on; grid on; grid minor;
plot(T,S(:,4),'b-','LineWidth',1.5);
plot(t_vec,x_sim(:,4),'r--','LineWidth',1.5);
ylim([-10 10]);
xlabel('time [sec]');
ylabel('Y velocity [km/s]');
xlim([0 P]);

saveas(fig,'ASEN5044_HW8_P2_sim.png','png');



function [ F ] = F_variant(X,Y)

mu = 398600;        % km^3/s^2
r0_nom = 6678;          % km
dt = 10;

F = expm(dt*[0, 1, 0, 0;
        mu*(2*X^2-Y^2)*(r0_nom)^(-5), 0, (3*mu*X*Y)*(r0_nom)^(-5), 0;
        0, 0, 0, 1;
        (3*mu*X*Y)*(r0_nom)^(-5), 0, -mu*(X^2-2*Y^2)*(r0_nom)^(-5), 0]);

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

