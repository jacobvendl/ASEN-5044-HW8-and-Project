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

%====== a ======%
Atil = [0, 1, 0, 0;
        2*mu/r0^3, 0, 0, 0;
        0, 0, 0, 1;
        0, 0, -mu/r0^3, 0];
    
Btil = [0, 0;
        1, 0;
        0, 0;
        0, 1];
    
%Ctil = [2*sqrt(r0-xs), 0, 2*sqrt(-ys), 0;
        


%====== b ======%




%====== c ======%
% One Full Orbit Period
% Simulate the Linearized DT dynamics
P = 2*pi*sqrt(r0^3/mu);

s0 = x0 + [1, 1, 0.1, 0.1]';
opts = odeset('RelTol',1e-12,'AbsTol',1e-12);
[T, S] = ode45(@(t,s)orbit_prop_func(t,s),[0 P],s0,opts);


fprintf('Plotting Satellite for 1 Orbit:\n');
fprintf("x0 = [6678, 0, 0, 7.7258]' + [1, 1, 0.1, 0.1]'");

fig = figure('visible','on');
set(fig,'Position',[100 100 900 600]);
sgtitle('Satellite State');

subplot(2,2,1); hold on; grid on; grid minor;
plot(T,S(:,1),'b-','LineWidth',1.5);
xlabel('time [sec]');
ylabel('X position [km]');
xlim([0 P]);

subplot(2,2,2); hold on; grid on; grid minor;
plot(T,S(:,2),'b-','LineWidth',1.5);
xlabel('time [sec]');
ylabel('Y position [km]');
xlim([0 P]);

subplot(2,2,3); hold on; grid on; grid minor;
plot(T,S(:,3),'b-','LineWidth',1.5);
xlabel('time [sec]');
ylabel('X velocity [km/s]');
xlim([0 P]);

subplot(2,2,4); hold on; grid on; grid minor;
plot(T,S(:,4),'b-','LineWidth',1.5);
xlabel('time [sec]');
ylabel('Y velocity [km/s]');
xlim([0 P]);

saveas(fig,'ASEN5044_HW8_P2_sim.png','png');








function [ ds ] = orbit_prop_func(t,s)

mu = 398600; 

x = s(1);
y = s(2);
r = sqrt(x^2+y^2);

xdot = s(3);
ydot = s(4);

xddot = -mu/r^3 * x;
yddot = -mu/r^3 * y;

ds = [xdot, ydot, xddot, yddot]';
end
