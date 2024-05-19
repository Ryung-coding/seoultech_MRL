%% ICCAS Simulation_Plant with Admiitance
clear 
clc
syms q1 q1_dot q1_ddot q2 q2_dot q2_ddot th1_d th2_d th1 th2 tau1 tau2 q1_pre q2_pre q1_dot_pre q2_dot_pre

% Plant Parameter
k1 = 15.05; % [N/mm]
k2 = 7.396; % [N/mm]
w_c = 30; % [rad/s]
L1 = 0.24215; % [m]
L2 = 0.1759;  % [m]
CoM1 = 0.13475; % [m]
CoM2 = 0.0562;  % [m]
m1 = 0.420;     % [kg]
m2 = 0.225;     % [kg]
J1 = m1*L1^2/12 + m1*(CoM1-L1/2)^2;  % [몰?루]
J2 = m2*L2^2/12 + m2*(L2/2-CoM2)^2;  % [몰?루]
g = 9.80665;

% Admittance Parameter
admit_m_x = 1;
admit_d_x = 1;
admit_k_x = 1;

admit_m_z = 1;
admit_d_z = 1;
admit_k_z = 1;

K_sp1 = 11.8;
K_sp2 = 5.7;

