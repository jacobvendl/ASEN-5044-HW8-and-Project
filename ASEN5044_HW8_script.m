%==========================================================================
% Jake Vendl | Jack Toland
% ASEN 5044
% Homework 8
% 12/3/2019
%==========================================================================
close all; clear all; clc

qw = 10;            % (m/s)^2
dt = 0.5;           % sec
OmegaA = 0.045;     % rad/s
OmegaB = -0.045;    % rad/s

FA = [1, sin(OmegaA*dt)/OmegaA, 0, -(1-cos(OmegaA*dt))/OmegaA;
    0, cos(OmegaA*dt), 0, -sin(OmegaA*dt);
    0, (1-cos(OmegaA*dt))/OmegaA, 1, sin(OmegaA*dt)/OmegaA;
    0, sin(OmegaA*dt), 0, cos(OmegaA*dt)];

FB = [1, sin(OmegaB*dt)/OmegaB, 0, -(1-cos(OmegaB*dt))/OmegaB;
    0, cos(OmegaB*dt), 0, -sin(OmegaB*dt);
    0, (1-cos(OmegaB*dt))/OmegaB, 1, sin(OmegaB*dt)/OmegaB;
    0, sin(OmegaB*dt), 0, cos(OmegaB*dt)];

AA = [0 1 0 0;
    0 0 0 -OmegaA;
    0 0 0 1;
    0 OmegaA 0 0];

AB = [0 1 0 0;
    0 0 0 -OmegaB;
    0 0 0 1;
    0 OmegaB 0 0];

%Problem 1

%====== a ======%
GammaA = [0 0; 1 0; 0 0; 0 1];
GammaB = [0 0; 1 0; 0 0; 0 1];

W = qw * [2 0.05; 0.05 0.5];

%find QA and QB using Van Loan's method
%first, find QA
ZA = dt.*[-AA GammaA*W*GammaA';
    zeros(4,4) AA'];
ezA = expm(ZA);
QA = ezA(5:8,5:8)' * ezA(1:4,5:8); %Q = (F')' * (inv(F)*Q)

ZB = dt.*[-AB GammaB*W*GammaB';
    zeros(4,4) AB'];
ezB = expm(ZB);
QB = ezB(5:8,5:8)' * ezB(1:4,5:8);

%====== b ======%
rng(100);

%b)i
H = [1 0 0 0; 0 0 1 0];
RA = [20 0.05; 0.05 20];
p=2;
T=200;

load('hw8problem1_data.mat');

Sv = chol(RA,'lower');
q = zeros(p,T); yA = zeros(p,T);
for k = 1:T
    qk = randn(2,1);
    vk = (Sv*qk);
    yA(:,k) = H*xasingle_truth(:,k) + vk;
end
figure; hold on; grid on; grid minor;
plot(1:1:40,yA(1,1:40),'-')
plot(1:1:40,yA(2,1:40),'-')
legend('yA_1','yA_2')
xlabel('time [s]'); ylabel('2D pseudo-measurements [m]')
title('(b.i) simulated measurements')
xticks([0 10 20 30 40])
xticklabels({'0','5','10','15','20'})

%b)ii
muA = [0, 85*cos(pi/4), 0, -85*sin(pi/4)]';
PA = 900*diag([10,2,10,2]);
% IMPLEMENT KALMAN HERE %
%define initial conditions
x_plus(:,1) = muA;

%handle the first time step outside of the loop, for indexing purposes
x_minus(:,1) = FA*muA;
P_minus(:,:,1) = FA*PA*FA' + QA;
P_plus(:,:,1) = PA;

%now iterate through the k's
for k=1:(T-1)
    x_minus(:,k+1) = FA*x_plus(:,k);
    P_minus(:,:,k+1) = FA*P_plus(:,:,k)*FA' + QA;
    K(:,:,k+1) = P_minus(:,:,k+1) * H' * inv(H*P_minus(:,:,k+1)*H'+RA);
    
    x_plus(:,k+1) = x_minus(:,k+1) + K(:,:,k+1) * (yA(:,k+1)-H*x_minus(:,k+1));
    P_plus(:,:,k+1) = (eye(4)-K(:,:,k+1)*H)*P_minus(:,:,k+1);
    
    sigma(1,k+1) = 2*sqrt(P_plus(1,1,k+1));
    sigma(2,k+1) = 2*sqrt(P_plus(2,2,k+1));
    sigma(3,k+1) = 2*sqrt(P_plus(3,3,k+1));
    sigma(4,k+1) = 2*sqrt(P_plus(4,4,k+1));
end

figure; hold on;
for i=1:4
    subplot(4,1,i); hold on; grid on; grid minor;
    plot(1:1:200,x_plus(i,:),'b-')
    plot(1:1:200,x_plus(i,:)+sigma(i,:),'k--')
    plot(1:1:200,x_plus(i,:)-sigma(i,:),'k--')
    xticks([0 20 40 60 80 100 120 140 160 180 200])
    xticklabels({'0','10','20','30','40','50','60','70','80','90','100'})
end
xlabel('time [s]');
subplot(4,1,1); ylabel('\xi [m]')
subplot(4,1,2); ylabel('\xiDot [m/s]')
subplot(4,1,3); ylabel('eta [m]')
subplot(4,1,4); ylabel('\etaDot [m/s]')
legend('component estimate','+/- 2\sigma')
suptitle('(b.ii) Kalman on aircraft A')


%====== c ======%
muB = [3200, 85*cos(pi/4), 3200, -85*sin(pi/4)]';
PB = 900*diag([11,4,11,4]);
RD = [10 0.15;
    0.15 10];

%c)i
FS = [FA zeros(4);
    zeros(4) FB];
AS = [AA zeros(4);
    zeros(4) AB];
WS = [W zeros(2);
    zeros(2) W];
GammaS = [GammaA zeros(4,2);
    zeros(4,2) GammaB];
ZS = dt*[-AS GammaS*WS*GammaS';
    zeros(8) AS'];
eZS = expm(ZS);
QS = eZS(9:16,9:16)' * eZS(1:8,9:16);

HS = [1 0 0 0 0 0 0 0;
    0 0 1 0 0 0 0 0;
    1 0 0 0 -1 0 0 0;
    0 0 1 0 0 0 -1 0];
RS = [RA zeros(2);
    zeros(2) RD];

%simulate measurements for yAprime and yD
SVA = chol(RA,'lower');
yAprime = zeros(p,T);
for k = 1:T
    qk = randn(2,1);
    vk = (SVA*qk);
    yAprime(:,k) = H*xadouble_truth(:,k) + vk;
end

SVD = chol(RD,'lower');
yD = zeros(p,T);
for k=1:T
    qk = randn(2,1);
    vk = (SVD*qk);
    yD(:,k) = [xadouble_truth(1,k); xadouble_truth(3,k)] - [xbdouble_truth(1,k); xbdouble_truth(3,k)] + vk;
end

yS = [yAprime; yD];

%now we are ready to get filtering
clear x_plus x_minus P_minus P_plus K

%define initial time step conditions prior to entering the loop
x_plus(:,1) = [muA; muB];
x_minus(:,1) = FS*x_plus(:,1);
P_minus(:,:,1) = FS*[PA zeros(4); zeros(4) PB]*FS' + QS;
P_plus(:,:,1) = [PA zeros(4); zeros(4) PB];

for k=1:(T-1)
    x_minus(:,k+1) = FS*x_plus(:,k);
    P_minus(:,:,k+1) = FS*P_plus(:,:,k)*FS' + QS;
    K(:,:,k+1) = P_minus(:,:,k+1) * HS' * inv(HS*P_minus(:,:,k+1)*HS'+RS);
    
    x_plus(:,k+1) = x_minus(:,k+1) + K(:,:,k+1) * (yS(:,k+1)-HS*x_minus(:,k+1));
    P_plus(:,:,k+1) = (eye(8)-K(:,:,k+1)*HS)*P_minus(:,:,k+1);
    
    sigma(1,k+1) = 2*sqrt(P_plus(1,1,k+1));
    sigma(2,k+1) = 2*sqrt(P_plus(2,2,k+1));
    sigma(3,k+1) = 2*sqrt(P_plus(3,3,k+1));
    sigma(4,k+1) = 2*sqrt(P_plus(4,4,k+1));
    sigma(5,k+1) = 2*sqrt(P_plus(5,5,k+1));
    sigma(6,k+1) = 2*sqrt(P_plus(6,6,k+1));
    sigma(7,k+1) = 2*sqrt(P_plus(7,7,k+1));
    sigma(8,k+1) = 2*sqrt(P_plus(8,8,k+1));
end

figure; hold on;
for i=1:4
    subplot(2,2,i); hold on; grid on; grid minor;
    plot(1:1:200, x_plus(2*i-1,:),'b-')
    plot(1:1:200, x_plus(2*i-1,:) + sigma(2*i-1,:),'k--')
    plot(1:1:200, x_plus(2*i-1,:) - sigma(2*i-1,:),'k--')
    xticks([0 20 40 60 80 100 120 140 160 180 200])
    xticklabels({'0','10','20','30','40','50','60','70','80','90','100'})
end
subplot(2,2,1); ylabel('\xi_A [m]'); legend('state component','+/- 2\sigma')
subplot(2,2,2); ylabel('\eta_A [m]'); xlabel('time [s]');
subplot(2,2,3); ylabel('\xi_B [m]'); xlabel('time [s]');
subplot(2,2,4); ylabel('\eta_B [m]'); xlabel('time [s');
suptitle('c)i position of aircraft A and B with time')

