%==========================================================================
% Jake Vendl | Jack Toland
% ASEN 5044
% Homework 8
% 12/3/2019
%==========================================================================

qw = 10;            % (m/s)^2
dt = 0.5;           % sec
OmegaA = 0.045;     % rad/s
OmegaB = -0.045;    % rad/s

F = [1, sin(Omega*dt)/Omega, 0, -(1-cos(Omega*dt))/Omega;
     0, cos(Omega*dt), 0, -sin(Omega*dt);
     0, (1-cos(Omega*dt))/Omega, 1, sin(Omega*dt)/Omega;
     0, sin(Omega*dt), 0, cos(Omega*dt)];


%% Problem 1

%====== a ======%
GammaA = [0 0; 1 0; 0 0; 0 1];
GammaB = [0 0; 1 0; 0 0; 0 1];

W = qw * [2 0.05; 0.05 0.5];


%====== b ======%
rng(100);

H = [1 0 0 0; 0 0 1 0];
RA = [20 0.05; 0.05 20];

load('hw8problem1data.mat');

% Sv = chol(R,'lower');
% q = zeros(p,T); y = zeros(p,T);
% for k = 1:T
%     q(:,k) = randn(p,1);
%     y(:,k) = x0 + Sv*q(:,k);
% end

muA = [0, 85*cos(pi/4), 0, -85*sin(pi/4)]';
PA = 900*diag([10,2,10,2]);

% IMPLEMENT KALMAN HERE %


%====== c ======%
muB = [3200, 85*cos(pi/4), 3200, -85*sin(pi/4)];
PB = 900*diag([11,4,11,4]);




%% Problem 2
% Orbit Determination %
%====== a ======%



%====== b ======%




%====== c ======%
% One Full Orbit Period
% Simulate the Linearized DT dynamics
