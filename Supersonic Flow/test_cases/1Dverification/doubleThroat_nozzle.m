clc;
clear;
close all;

% Important Constants
gam = 1.4; % <- might not be the right constant?
M_air = 28.9647/1000; % kg/mol

%% Declare Inputs

% Inlet Conditions
P_i = 50e3; % pascals
T_i = 368.16; %K
u_i = 5.42; %m/s

% Exit Conditions
P_e = 15e3; % Pascals

% Nozzle Geometry functions
aftThrArea = @(x) pi.*(-sin(10.*pi.*x)./250 + x./100 + 1/100).^2;
foreThrArea = @(x) pi.*(-cos(10.*pi.*(x-0.05))/50 + 0.0265).^2;

%% Get Inlet Mach Number and Other Properties

AR_inlet = foreThrArea(-0.05)/aftThrArea(0.05);
% function comes from Elements of Propulsion by Mattingly
func = @(M) (1./M).*((2/(gam+1)).*(1+0.5.*(gam-1).*M.^2)).^((gam+1)/(2*(gam-1))) - AR_inlet;
jFunc = @(M) (2.*(M.^2-1).*(((gam-1).*M.^2+1)./(gam+1)).^(1/(gam-1) - 0.5))./((gam+1).*M.^2);
M_init = 0.05;
tol = 1e-7;

% Iterate for Mach Number using Newton's method
M_i = newtonSys(func, jFunc, M_init, tol, tol, 1e5, 1);

% Check if Gas is Air
% a_i = u_i/M_i;
R_gas = 287.058;
% R_gas = (a_i^2)./(gam*T_i); % (m^2/s^2)/K -> J/(kg*K)
a_i = sqrt(gam*R_gas*T_i);
H_i = (a_i^2)/(gam-1) + 0.5.*(M_i*a_i)^2; % m^2 / s^2

% Get Density
rho_i = gam.*P_i./(a_i^2); % kg/m^3

% Mass flow
m_i = rho_i*foreThrArea(-0.05)*(M_i*a_i);

P_total = P_i*(1+0.5*(gam-1)*M_i^2)^(gam/(gam-1));

%% Calculate Throat Conditions

a_T = sqrt(H_i*2*(gam-1)/(gam+1));
rho_T = rho_i*u_i*AR_inlet / a_T;
P_T = (a_T^2)*rho_T/gam;
T_T = P_T./(rho_T*R_gas);
m_T = rho_T*aftThrArea(0.05)*a_T;

P_T_total = P_T*(1+0.5*(gam-1))^(gam/(gam-1));

%% Calculate Exit Conditions

% Calculate Exit mach Number
AR_exit = aftThrArea(0.35)/aftThrArea(0.05);
func = @(M) (1./M).*((2/(gam+1)).*(1+0.5.*(gam-1).*M.^2)).^((gam+1)/(2*(gam-1))) - AR_exit;
jFunc = @(M) (2.*(M.^2-1).*(((gam-1).*M.^2+1)./(gam+1)).^(1/(gam-1) - 0.5))./((gam+1).*M.^2);
M_e = newtonSys(func, jFunc, 1.6, tol, tol, 1e5, 1);

% Calculate Using Exit Pressures
p_e_ratio = P_e/(rho_T*a_T^2);

AA = 0.5;
CC = -3;
BB = (gam/(gam-1))*p_e_ratio*AR_exit;

u_1 = (-BB + sqrt(BB^2 - 4*AA*CC))./(2*AA);
u_2 = (-BB - sqrt(BB^2 - 4*AA*CC))./(2*AA);

%% Print Out Results?

fprintf('Inlet Conditions: \n');
fprintf('Inlet Pressure: %0.2f Pa\n', P_i);
fprintf('Inlet Temperature: %0.2f K\n', T_i);
fprintf('Inlet Density: %0.2f kg/m^3\n', rho_i);
fprintf('Inlet Speed: %0.2f m/s\n', M_i*a_i);
fprintf('Inlet Sound Speed: %0.2f m/s\n', a_i);
fprintf('Inlet Mach Number: %0.3f\n', M_i);
fprintf('Inlet Mass Flow: %0.3f kg/s\n', m_i);
fprintf('Gas Constant: %0.3f J/(mol * K)\n', R_gas*M_air); % <- not air???


fprintf('\nThroat Conditions: \n');
fprintf('Throat Pressure: %0.2f Pa\n', P_T);
fprintf('Throat Sound Speed: %0.2f m/s\n', a_T);
fprintf('Throat Density: %0.3f kg/m^3\n', rho_T);
fprintf('Throat Temperature: %0.3f K\n', T_T);
fprintf('Throat Mass Flow: %0.3f kg/s\n', m_T);

save('doubleThroat', 'p_e_ratio', 'aftThrArea', 'u_1');