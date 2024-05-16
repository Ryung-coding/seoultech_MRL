clear all
clc

%envirment setting
Ts = 0.01; %s
sampling_num = 500;
N = 10; % Prediction horizon


xr_q = [deg2rad(20); 0; deg2rad(20); 0];
xr = [xr_q; 0;0;0;0];
x0_q = [0; 0; 0; 0];
x0 = [x0_q; 0;0;0;0];
x0_initial=x0;

q1=zeros(1,sampling_num);
q1_dot=zeros(1,sampling_num);
q2=zeros(1,sampling_num);
q2_dot=zeros(1,sampling_num);

theta1=zeros(1,sampling_num);
theta1_dot=zeros(1,sampling_num);
theta2=zeros(1,sampling_num);
theta2_dot=zeros(1,sampling_num);

theta1_des=zeros(1,sampling_num);
theta2_des=zeros(1,sampling_num);

X = [q1; q1_dot; q2; q2_dot; theta1; theta1_dot; theta2; theta2_dot];
X(:,1) = x0;
U = [theta1_des; theta2_des];


% SEA state space
J1 = 0.420*0.24215^2/12+0.420*0.13475^2;
J2 = 0.225*0.1759^2/12+0.225*0.0562^2;
k1 = 15.75;
k2 = 6.75;
c1 = 0.01;
c2 = 0.01;

A_SEA = [     0               1           0          0;
         -(k1+k2)/J1    -(c1+c2)/J1     k2/J1      c2/J1;
              0               0           0          1;
            k2/J2           c2/J2      -k2/J2     -c2/J2];

B_SEA = [  0      0        0         0;
         k1/J1  c1/J1   -k2/J1    -c2/J1;
           0      0        0         0;
           0      0      k2/J2    c2/J2];

C_SEA = eye(size(A_SEA));

D_SEA = 0;
sys_SEA = ss(A_SEA, B_SEA, C_SEA, D_SEA);
sysd_SEA = c2d(sys_SEA, Ts, 'zoh');
Ad_SEA=sysd_SEA.A;
Bd_SEA=sysd_SEA.B;

%% Bessel Filter state space
w0_1=7;
w0_2=7;

sys_Bessel=tf({[0 0 3] [0 0 3]},{[1/w0_1^2 3/w0_1 3] [1/w0_2^2 3/w0_2 3]})
%bode(sys_Bessel)

% x=[theta1 theta1_dot theta2 theta2_dot]
% u=[theta1_des theta2_des]
A_motor=[     0        1        0           0;
         -3*w0_1^2  -3*w0_1     0           0;
              0        0        0           1;
              0        0    -3*w0_2^2    -3*w0_2];

B_motor=[    0       0;
          3*w0_1^2   0;
             0       0;
             0    3*w0_2^2];

C_motor = eye(size(A_motor));

D_motor = 0;
sys_motor = ss(A_motor, B_motor, C_motor, D_motor);
sysd_motor = c2d(sys_motor, Ts, 'zoh');
Ad_motor=sysd_motor.A;
Bd_motor=sysd_motor.B;

buf_x=zeros(4,500);
buf_x(:,1)=[-1 0 1 0];
for k=1:1:499
    buf_x(:,k+1)=Ad_motor*buf_x(:,k);
end
plot(buf_x(3,:))



%% butterworth filter state space
zeta = 0.5;
w0_1=5;
w0_2=5;

sys_butterworth = tf({[0 0 w0_1^2] [0 0 w0_2^2]},{[1 2*zeta*w0_1 w0_1^2] [1 2*zeta*w0_2 w0_2^2]})
%bode(sys_butterworth)

% x=[theta1 theta1_dot theta2 theta2_dot]
% u=[theta1_des theta2_des]
A_motor=[     0          1          0             0;
           -w0_1^2  -2*zeta*w0_1    0             0;
              0          0          0             1;
              0          0       -w0_2^2    -2*zeta*w0_2];

B_motor=[    0       0;
            w0_1^2   0;
             0       0;
             0     w0_2^2];

C_motor = eye(size(A_motor));

D_motor = 0;
sys_motor = ss(A_motor, B_motor, C_motor, D_motor);
sysd_motor = c2d(sys_motor, Ts, 'zoh');
Ad_motor=sysd_motor.A;
Bd_motor=sysd_motor.B;

buf_x=zeros(4,500);
buf_x(:,1)=[-1 0 1 0];
for k=1:1:499
    buf_x(:,k+1)=Ad_motor*buf_x(:,k);
end
plot(buf_x(3,:))


%% SOLVER setting

umin = [-(x0(1)+deg2rad(100)) ; -(x0(3)+deg2rad(100))];
umax = [ (x0(1)+deg2rad(100)) ;  (x0(3)+deg2rad(100))];
xmin = [-deg2rad(100) ; -deg2rad(100) ; -deg2rad(100) ; -deg2rad(100) ; umin ; umin];
xmax = [ deg2rad(100) ;  deg2rad(100) ;  deg2rad(100) ;  deg2rad(100) ; umax ; umax];

% Objective function
Q = diag([100 1 10 1 0 0 0 0]);       %[q1 q1_dot q2 q2_dot theta1 theta1_dot theta2 theta2_dot]
QN = Q;                              %[q1 q1_dot q2 q2_dot theta1 theta1_dot theta2 theta2_dot]

R_theta_des = diag([0.01 0.01]);         %[theta1_des theta2_des]
R = R_theta_des;                     %[theta1_des theta2_des]

Ad_model=[       Ad_SEA                  Bd_SEA;
                zeros(size(Ad_motor))   Ad_motor];

Bd_model=[       zeros(size(Bd_SEA,1),2); 
                Bd_motor              ];

[nx, nu] = size(Bd_model);

%%
P = blkdiag( kron(speye(N), Q), QN, kron(speye(N), R) );
% - linear objective
q = [repmat(-Q*xr, N, 1); -QN*xr; zeros(N*nu, 1)];
% - linear dynamics
Ax = kron(speye(N+1), -speye(nx)) + kron(sparse(diag(ones(N, 1), -1)), Ad_model);
Bu = kron([sparse(1, N); speye(N)], Bd_model);
Aeq = [Ax, Bu];
leq = [-x0; zeros(N*nx, 1)];
ueq = leq;
% - input and state constraints
Aineq = speye((N+1)*nx + N*nu);
lineq = [repmat(xmin, N+1, 1); repmat(umin, N, 1)];
uineq = [repmat(xmax, N+1, 1); repmat(umax, N, 1)];
% - OSQP constraints
A = [Aeq; Aineq];
l = [leq; lineq];
u = [ueq; uineq];

% Create an OSQP object
prob = osqp;

% Setup workspace
prob.setup(P, q, A, l, u);

% Simulate in closed loop
for i = 1 : sampling_num
    % Solve
    res = prob.solve();

    % Check solver status
    if ~strcmp(res.info.status, 'solved')
        error('OSQP did not solve the problem!')
    end

    % Apply first control input to the plant
    ctrl = res.x((N+1)*nx+1:(N+1)*nx+nu);

    x0 = Ad_model*x0 + Bd_model*(ctrl);

    X(:,i+1) = x0;
    U(:,i+1) = ctrl;

    
    % Update
    umin = [-(x0(1)+deg2rad(100)) ; -(x0(3)+deg2rad(100))];
    umax = [ (x0(1)+deg2rad(100)) ;  (x0(3)+deg2rad(100))];
    l(1:nx) = -x0;
    %l((N+1)*nx+1:end) = [repmat(xmin, N+1, 1); repmat(umin, N, 1)];
    u(1:nx) = -x0;
    %u((N+1)*nx+1:end) = [repmat(xmax, N+1, 1); repmat(umax, N, 1)];
    prob.update('l', l, 'u', u);


end

% Assuming data preparation for Xr and X0
Xr = xr*ones(1, sampling_num); % Replicating Xr sampling_num times
X0 = x0_initial*ones(1, sampling_num); % Replicating X0 sampling_num times

set(gcf, 'Color', 'w'); % Setting the background to white
% Drawing each subplot
subplot(5,2,1)
plot(X(1,:), 'LineWidth', 2)
hold on
plot(1, X0(1,1), 'bo', 'MarkerSize', 6, 'LineWidth', 1);
plot(Xr(1,:), '--r', 'LineWidth', 3)
hold off
ylim([-1 1])
xlabel('step', 'FontName', 'Times New Roman', 'FontSize', 12, 'FontWeight', 'bold')
ylabel('q1', 'FontName', 'Times New Roman', 'FontSize', 12, 'FontWeight', 'bold')

subplot(5,2,2)
plot(X(2,:), 'LineWidth', 2)
hold on
plot(1, X0(2,1), 'bo', 'MarkerSize', 6, 'LineWidth', 1);
plot(Xr(2,:), '--r', 'LineWidth', 2)
hold off
ylim([-1 1])
xlabel('step', 'FontName', 'Times New Roman', 'FontSize', 12, 'FontWeight', 'bold')
ylabel('q1_dot', 'FontName', 'Times New Roman', 'FontSize', 12, 'FontWeight', 'bold')

subplot(5,2,3)
plot(X(3,:), 'LineWidth', 2)
hold on
plot(1, X0(3,1), 'bo', 'MarkerSize', 6, 'LineWidth', 1);
plot(Xr(3,:), '--r', 'LineWidth', 2)
hold off
ylim([-1 1])
xlabel('step', 'FontName', 'Times New Roman', 'FontSize', 12, 'FontWeight', 'bold')
ylabel('q2', 'FontName', 'Times New Roman', 'FontSize', 12, 'FontWeight', 'bold')

subplot(5,2,4)
plot(X(4,:), 'LineWidth', 2)
hold on
plot(1, X0(4,1), 'bo', 'MarkerSize', 6, 'LineWidth', 1);
plot(Xr(4,:), '--r', 'LineWidth', 2)
hold off
ylim([-1 1])
xlabel('step', 'FontName', 'Times New Roman', 'FontSize', 12, 'FontWeight', 'bold')
ylabel('q2_dot', 'FontName', 'Times New Roman', 'FontSize', 12, 'FontWeight', 'bold')

subplot(5,2,5)
plot(X(5,:), 'LineWidth', 2)
hold on
plot(1, X0(5,1), 'bo', 'MarkerSize', 6, 'LineWidth', 1);
plot(Xr(5,:), '--r', 'LineWidth', 3)
hold off
ylim([-1 1])
xlabel('step', 'FontName', 'Times New Roman', 'FontSize', 12, 'FontWeight', 'bold')
ylabel('theta1', 'FontName', 'Times New Roman', 'FontSize', 12, 'FontWeight', 'bold')

subplot(5,2,6)
plot(X(6,:), 'LineWidth', 2)
hold on
plot(1, X0(6,1), 'bo', 'MarkerSize', 6, 'LineWidth', 1);
plot(Xr(6,:), '--r', 'LineWidth', 3)
hold off
ylim([-1 1])
xlabel('step', 'FontName', 'Times New Roman', 'FontSize', 12, 'FontWeight', 'bold')
ylabel('theta1_dot', 'FontName', 'Times New Roman', 'FontSize', 12, 'FontWeight', 'bold')

subplot(5,2,7)
plot(X(7,:), 'LineWidth', 2)
hold on
plot(1, X0(7,1), 'bo', 'MarkerSize', 6, 'LineWidth', 1);
plot(Xr(7,:), '--r', 'LineWidth', 3)
hold off
ylim([-1 1])
xlabel('step', 'FontName', 'Times New Roman', 'FontSize', 12, 'FontWeight', 'bold')
ylabel('theta2', 'FontName', 'Times New Roman', 'FontSize', 12, 'FontWeight', 'bold')

subplot(5,2,8)
plot(X(8,:), 'LineWidth', 2)
hold on
plot(1, X0(8,1), 'bo', 'MarkerSize', 6, 'LineWidth', 1);
plot(Xr(8,:), '--r', 'LineWidth', 3)
hold off
ylim([-1 1])
xlabel('step', 'FontName', 'Times New Roman', 'FontSize', 12, 'FontWeight', 'bold')
ylabel('theta2_dot', 'FontName', 'Times New Roman', 'FontSize', 12, 'FontWeight', 'bold')

subplot(5,2,9)
plot(U(1,:), 'LineWidth', 2)
ylim([-1 1])
xlabel('step', 'FontName', 'Times New Roman', 'FontSize', 12, 'FontWeight', 'bold')
ylabel('theta1_des', 'FontName', 'Times New Roman', 'FontSize', 12, 'FontWeight', 'bold')

subplot(5,2,10)
plot(U(2,:), 'LineWidth', 2)
ylim([-1 1])
xlabel('step', 'FontName', 'Times New Roman', 'FontSize', 12, 'FontWeight', 'bold')
ylabel('theta2_des', 'FontName', 'Times New Roman', 'FontSize', 12, 'FontWeight', 'bold')





%%

Wc = ctrb(Ad_model,Bd_model);
rankWc = rank(Wc)

if rankWc == size(A,1)
    disp('System is controllable.');
else
    disp('System is not controllable.');
end

Wo = obsv(Ad_model,eye(size(Ad_model)));
rankWo = rank(Wo)

if rankWo == size(A,1)
    disp('System is observable.');
else
    disp('System is not observable.');
end
