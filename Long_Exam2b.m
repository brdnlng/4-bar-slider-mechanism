% MAE 364-01
% Spring 2021
% Exam 2b
%
% From: Braeden Long
% Date: 04/05/2021

%% Defining variables and vectors

clc; clear; close all
global r1 r2 theta2 omega2 r3 velocity3 acceleration3 r4 theta4 omega4 alpha4
r1 = 10;                % length of ground, cm
r2 = 4;                 % length of driver, cm
r4 = 4;                 % length of coupler, cm
omega2 = 15*(2*pi/60);  % constant angular velocity of crank (counter-clockwise), rad/s
D2R = pi/180;            % convert degrees to radians
R2D = 1/D2R;              % convert radians to degrees
inc = 0:2:360;          % array of degree values, increments of 2
Vtheta2 = inc*D2R;       % vector of theta2 values in radians 

% predefined vectors for values needing to be determined for each increment
Vtheta4 = zeros(1,length(Vtheta2));
Vr3 = Vtheta4;
Vomega4 = Vtheta4;
Vvelocity3 = Vtheta4;
Valpha4 = Vtheta4;
Vacceleration3 = Vtheta4;

% initial guesses when theta2 = 0
theta4 = 135*D2R; r3 = 10;
xInital = [theta4, r3];     % positions
vInital = [0 0];            % velocities
aInital = [0 0];            % accelerations

%% Loop to determine unknowns

opt = optimset('Display', 'off');
for i = 1:length(Vtheta2) % loop for every value of theta2
  theta2 = Vtheta2(i);
% Solving for position using @PosEq
  [Xtemp, ~] = fsolve(@PosEq,xInital,opt);
  theta4 = Xtemp(1);    % solution for theta4 for current theta2 value
  r3 = Xtemp(2);        % solution for r3 for current theta2 value
  Vtheta4(i) = theta4;  % storing theta4 solution
  Vr3(i) = r3;          % storing r3 solution
  xInital = [theta4,r3];

% Solving for angular velocities using @VelEq
  [Vtemp, ~] = fsolve(@VelEq,vInital,opt);
  omega4 = Vtemp(1);            % solution for omega4 for current theta2 value
  velocity3 = Vtemp(2);         % solution for velocity3 for current theta2 value
  Vomega4(i)= omega4;           % storing omega4 solution
  Vvelocity3(i) = velocity3;    % storing velocity3 solution
  vInital = [omega4, velocity3];
 
% Solving for angular accelerations using @AccEq
  [Atemp, ~] = fsolve(@AccEq,aInital,opt);
  alpha4 = Atemp(1);                    % solution for alpha4 for current theta2 value
  acceleration3 = Atemp(2);             % solution for acceleration3 for current theta2 value
  Valpha4(i)= alpha4;                   % storing alpha4 solution
  Vacceleration3(i) = acceleration3;    % storing acceleration3 solution
  aInital = [alpha4, acceleration3];
end

%% Plots

% plot length, velocity, and acceleration for link 3 versus theta2
figure(1)
subplot(3,1,1)
  plot(R2D*Vtheta2,Vr3, '-r');
  xlim([0 360])
  ylabel('r_3 (cm)');
  set(gca,'xtick',0:30:360)
subplot(3,1,2)
  plot(R2D*Vtheta2,Vvelocity3, '-g');
  xlim([0 360])
  ylabel('v_3 (cm/s)');
  set(gca,'xtick',0:30:360)
subplot(3,1,3)
  plot(R2D*Vtheta2,Vacceleration3, '-b');
  xlim([0 360])
  ylabel('a_3 (cm/s^2)');
  xlabel('\theta_2');
  set(gca,'xtick',0:30:360)

% plot angular position, velocity, acceleration of link 4 versus theta2
figure(2)
subplot(3,1,1)
  plot(R2D*Vtheta2,R2D*Vtheta4, '-c');
  xlim([0 360])
  ylabel('\theta_4');
  set(gca,'xtick',0:30:360)
subplot(3,1,2)
  plot(R2D*Vtheta2,R2D*Vomega4, '-m');
  xlim([0 360])
  ylabel('\omega_4');
  set(gca,'xtick',0:30:360)
subplot(3,1,3)
  plot(R2D*Vtheta2,R2D*Valpha4, '-y');
  xlim([0 360])
  ylabel('\alpha_4');
  xlabel('\theta_2');
  set(gca,'xtick',0:30:360)
  
%% Animation

% x,y coordinates for origin
Ox = zeros(1,length(Vtheta2));
Oy = zeros(1,length(Vtheta2));

% x,y coordinates for r1
r1x = zeros(1,length(Vtheta2));
r1y = ones(1,length(Vtheta2))*10;

% x,y coordinates for r2
r2x = real(r2*exp(1i*Vtheta2)) + r1x;
r2y = imag(r2*exp(1i*Vtheta2)) + r1y;

% x,y coordinates for r4
r4x = real(r4*exp(1i*Vtheta4));
r4y = imag(r4*exp(1i*Vtheta4));

% x,y coordinates for r5
r5x = zeros(1,length(Vtheta2));
r5y = zeros(1,length(Vtheta2));
for i = 1:length(Vtheta2)
    r5x(i) = real((16-Vr3(i))*exp(1i*(Vtheta4(i)-pi/2))) + r2x(i);
    r5y(i) = imag((16-Vr3(i))*exp(1i*(Vtheta4(i)-pi/2))) + r2y(i);
end

% x,y coordinates for r6
r6x = real(8*exp(1i*Vtheta4)) + r5x;
r6y = imag(8*exp(1i*Vtheta4)) + r5y;

% Mechanism Animation
figure
for i = 1:length(Vtheta2)
    bar1x = [Ox(i) r1x(i)];     % Coordinates Link 1
    bar1y = [Oy(i) r1y(i)];
    bar2x = [r1x(i) r2x(i)];    % Coordinates Link 2
    bar2y = [r1y(i) r2y(i)];
    bar4x = [Ox(i) r4x(i)];     % Coordinates Link 4
    bar4y = [Oy(i) r4y(i)];
    bar3x = [r4x(i) r2x(i)];    % Coordinates Link 3
    bar3y = [r4y(i) r2y(i)];
    bar5x = [r2x(i) r5x(i)];    % Coordinates Link 5
    bar5y = [r2y(i) r5y(i)];
    bar6x = [r5x(i) r6x(i)];    % Coordinates Link 6
    bar6y = [r5y(i) r6y(i)];    
    
    plot(bar1x,bar1y,bar2x,bar2y,bar3x,bar3y,bar4x,bar4y,bar5x,bar5y,bar6x,bar6y);
    axis([-15 15 -5 25]);
    M(i)=getframe; % For assembling movie frames
end

%% Functions

function F = PosEq(X)
% function to solve for position equations
global theta2 theta4 r1 r2 r3 r4
theta4 = X(1); r3 = X(2); % locally rename solution variables
% complex vector equation for position
f = 1i*r1 + r2*exp(1i*theta2) - r3*exp(1i*(theta4-pi/2)) - r4*exp(1i*theta4);
F = [real(f); imag(f)];
end
 
function F = VelEq(X)
% function to solve for velocity equations  
global theta2 r2 omega2 r3 r4 omega4 theta4 velocity3
omega4 = X(1); velocity3 = X(2); % locally rename solution variables
% complex vector equation for velocity
f = 1i*r2*omega2*exp(1i*theta2) - 1i*r4*omega4*exp(1i*theta4) - ...
    1i*r3*omega4*exp(1i*omega4-pi/2);
F = [real(f); imag(f)];
end

function F = AccEq(X)
% function to solve for acceleration equations  
global theta2 r2 r3 r4 theta4 omega2 omega4 alpha4 acceleration3 velocity3
alpha4 = X(1); acceleration3 = X(2);  % Locally rename solution variables
% complex vector equation for accleration
f = 1i*r2*exp(1i*theta2)*(1i*omega2^2) - 1i*r4*exp(1i*theta4)*(alpha4 + 1i*omega4^2) - ...
    exp(1i*(theta4-pi/2))*(acceleration3 + 1i*2*velocity3*omega4 + 1i*r3*alpha4 + 1i*r3*omega4^2);
F = [real(f); imag(f)];
end