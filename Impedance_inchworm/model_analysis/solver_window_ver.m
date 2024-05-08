
%envirment setting
Ts = 0.01; %s
tau_stall = 2.5; %Nm
r = 0.03; %m
sampling_num = 500;
N = 10; % Prediction horizon
x0 = [1; 0.7; -0.5; -1; 0; 0]; %initial state
x0_initial=x0;
%model parameter setting 
J1 = 0.420*0.24215^2/3;
J2 = 0.225*0.1759^2/3;

k1 = 2*9.8;
k2 = 2*5.9;

c1 = 0.0001;
c2 = 0.0001;

%python_ver
J1 = 1;
J2 = 1;

k1 = 0.1;
k2 = 0.1;

c1 = 0.001;
c2 = 0.001;
x0 = [0; 0; 0; 0; 0; 0]; %initial state
x0_initial=x0;


%control value
xr = [1; 0; 1; 0; 0; 0]; %rad rad/2 rad rad/2 Nm Nm

q1=zeros(1,sampling_num);
q1_dot=zeros(1,sampling_num);
q2=zeros(1,sampling_num);
q2_dot=zeros(1,sampling_num);
tau1=zeros(1,sampling_num);
tau2=zeros(1,sampling_num);
theta1=zeros(1,sampling_num);
theta2=zeros(1,sampling_num);

X = [q1; q1_dot; q2; q2_dot; tau1; tau2];
X(:,1) = x0;

U = [theta1; theta2];

A = [  0       1       0       0    0    0;
    -k1/J1  -c1/J1   k2/J1     0   1/J1  0;
       0       0       0       1    0    0;
       0       0    -k2/J2  -c2/J2  0   1/J2;
       0       0       0       0    0    0;
       0       0       0       0    0    0];

B = [  0        0;
     k1/J1   -k2/J1;
       0        0;
       0      k2/J2;
       0        0;
       0        0];

C = eye(size(A));

D = 0;

sys = ss(A, B, C, D);
sysd = c2d(sys, Ts, 'zoh');
Ad=sysd.A;
Bd=sysd.B;
[nx, nu] = size(Bd);

% Discrete time model of a quadcopter

% Constraints
umin = [-tau_stall/(r*k1); -tau_stall/(r*k2)];
umax = [ tau_stall/(r*k1);  tau_stall/(r*k2)];
xmin = [-pi/2; -pi/2; -pi/2; -pi/2; -Inf; -Inf];
xmax = [ pi/2;  pi/2;  pi/2;  pi/2;  Inf;  Inf];
% Constraints_py ver
d = 25.5*0.1;
umin = [-d*(100.*3.141592/180.); -d*(100.*3.141592/180.)];
umax = [ d*(100.*3.141592/180.);  d*(100.*3.141592/180.)];
xmin = [-(10.*3.141592/180.); -(90.*3.141592/180.); -(10.*3.141592/180.); -(90.*3.141592/180.); -Inf; -Inf];
xmax = [(10.*3.141592/180.); (90.*3.141592/180.); (10.*3.141592/180.); (90.*3.141592/180.); Inf; Inf];




% Objective function
Q = diag([10 10 10 10 0 0]); %position control version
QN = Q;
R = 0.05*eye(2);

% Cast MPC problem to a QP: x = (x(0),x(1),...,x(N),u(0),...,u(N-1))
% - quadratic objective
P = blkdiag( kron(speye(N), Q), QN, kron(speye(N), R) );
% - linear objective
q = [repmat(-Q*xr, N, 1); -QN*xr; zeros(N*nu, 1)];
% - linear dynamics
Ax = kron(speye(N+1), -speye(nx)) + kron(sparse(diag(ones(N, 1), -1)), Ad);
Bu = kron([sparse(1, N); speye(N)], Bd);
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

%% Simulate in closed loop
for i = 1 : sampling_num
    % Solve
    res = prob.solve();

    % Check solver status
    if ~strcmp(res.info.status, 'solved')
        error('OSQP did not solve the problem!')
    end

    % Apply first control input to the plant
    ctrl = res.x((N+1)*nx+1:(N+1)*nx+nu);

    x0 = Ad*x0 + Bd*(ctrl);

    X(:,i+1) = x0;
    U(:,i+1) = ctrl;


    % Update initial state
    l(1:nx) = -x0;
    u(1:nx) = -x0;
    prob.update('l', l, 'u', u);
    
    
end

%% Assuming data preparation for Xr and X0
Xr = xr*ones(1, sampling_num); % Replicating Xr sampling_num times
X0 = x0_initial*ones(1, sampling_num); % Replicating X0 sampling_num times

set(gcf, 'Color', 'w'); % Setting the background to white

% Drawing each subplot
subplot(4,2,1)
plot(X(1,:), 'LineWidth', 2)
hold on
plot(1, X0(1,1), 'bo', 'MarkerSize', 6, 'LineWidth', 1);
plot(Xr(1,:), '--r', 'LineWidth', 3)
hold off
ylim([-1 1])
xlabel('step', 'FontName', 'Times New Roman', 'FontSize', 12, 'FontWeight', 'bold')
ylabel('q1', 'FontName', 'Times New Roman', 'FontSize', 12, 'FontWeight', 'bold')

subplot(4,2,2)
plot(X(2,:), 'LineWidth', 2)
hold on
plot(1, X0(2,1), 'bo', 'MarkerSize', 6, 'LineWidth', 1);
plot(Xr(2,:), '--r', 'LineWidth', 2)
hold off
ylim([-1 1])
xlabel('step', 'FontName', 'Times New Roman', 'FontSize', 12, 'FontWeight', 'bold')
ylabel('q1_dot', 'FontName', 'Times New Roman', 'FontSize', 12, 'FontWeight', 'bold')

subplot(4,2,3)
plot(X(3,:), 'LineWidth', 2)
hold on
plot(1, X0(3,1), 'bo', 'MarkerSize', 6, 'LineWidth', 1);
plot(Xr(3,:), '--r', 'LineWidth', 2)
hold off
ylim([-1 1])
xlabel('step', 'FontName', 'Times New Roman', 'FontSize', 12, 'FontWeight', 'bold')
ylabel('q2', 'FontName', 'Times New Roman', 'FontSize', 12, 'FontWeight', 'bold')

subplot(4,2,4)
plot(X(4,:), 'LineWidth', 2)
hold on
plot(1, X0(4,1), 'bo', 'MarkerSize', 6, 'LineWidth', 1);
plot(Xr(4,:), '--r', 'LineWidth', 2)
hold off
ylim([-1 1])
xlabel('step', 'FontName', 'Times New Roman', 'FontSize', 12, 'FontWeight', 'bold')
ylabel('q2_dot', 'FontName', 'Times New Roman', 'FontSize', 12, 'FontWeight', 'bold')

subplot(4,2,5)
plot(X(5,:), 'LineWidth', 2)
hold on
plot(1, X0(5,1), 'bo', 'MarkerSize', 6, 'LineWidth', 1);
hold off
ylim([-1 1])
xlabel('step', 'FontName', 'Times New Roman', 'FontSize', 12, 'FontWeight', 'bold')
ylabel('tau1', 'FontName', 'Times New Roman', 'FontSize', 12, 'FontWeight', 'bold')

subplot(4,2,6)
plot(X(6,:), 'LineWidth', 2)
hold on
plot(1, X0(6,1), 'bo', 'MarkerSize', 6, 'LineWidth', 1);
hold off
ylim([-1 1])
xlabel('step', 'FontName', 'Times New Roman', 'FontSize', 12, 'FontWeight', 'bold')
ylabel('tau2', 'FontName', 'Times New Roman', 'FontSize', 12, 'FontWeight', 'bold')

subplot(4,2,7)
plot(U(1,:), 'LineWidth', 2)
ylim([-1 1])
xlabel('step', 'FontName', 'Times New Roman', 'FontSize', 12, 'FontWeight', 'bold')
ylabel('theta1', 'FontName', 'Times New Roman', 'FontSize', 12, 'FontWeight', 'bold')

subplot(4,2,8)
plot(U(2,:), 'LineWidth', 2)
ylim([-1 1])
xlabel('step', 'FontName', 'Times New Roman', 'FontSize', 12, 'FontWeight', 'bold')
ylabel('theta2', 'FontName', 'Times New Roman', 'FontSize', 12, 'FontWeight', 'bold')




%% example

Ad = [1       0       0   0   0   0   0.1     0       0    0       0       0;
      0       1       0   0   0   0   0       0.1     0    0       0       0;
      0       0       1   0   0   0   0       0       0.1  0       0       0;
      0.0488  0       0   1   0   0   0.0016  0       0    0.0992  0       0;
      0      -0.0488  0   0   1   0   0      -0.0016  0    0       0.0992  0;
      0       0       0   0   0   1   0       0       0    0       0       0.0992;
      0       0       0   0   0   0   1       0       0    0       0       0;
      0       0       0   0   0   0   0       1       0    0       0       0;
      0       0       0   0   0   0   0       0       1    0       0       0;
      0.9734  0       0   0   0   0   0.0488  0       0    0.9846  0       0;
      0      -0.9734  0   0   0   0   0      -0.0488  0    0       0.9846  0;
      0       0       0   0   0   0   0       0       0    0       0       0.9846];
Bd = [0      -0.0726  0       0.0726;
     -0.0726  0       0.0726  0;
     -0.0152  0.0152 -0.0152  0.0152;
      0      -0.0006 -0.0000  0.0006;
      0.0006  0      -0.0006  0;
      0.0106  0.0106  0.0106  0.0106;
      0      -1.4512  0       1.4512;
     -1.4512  0       1.4512  0;
     -0.3049  0.3049 -0.3049  0.3049;
      0      -0.0236  0       0.0236;
      0.0236  0      -0.0236  0;
      0.2107  0.2107  0.2107  0.2107]; 
[nx, nu] = size(Bd);

% Constraints
u0 = 10.5916;
umin = [9.6; 9.6; 9.6; 9.6] - u0;
umax = [13; 13; 13; 13] - u0;
xmin = [-pi/6; -pi/6; -Inf; -Inf; -Inf; -1; -Inf(6,1)];
xmax = [ pi/6;  pi/6;  Inf;  Inf;  Inf; Inf; Inf(6,1)];

% Objective function
Q = diag([0 0 10 10 10 10 0 0 0 5 5 5]);
QN = Q;
R = 0.1*eye(4);

% Initial and reference states
x0 = zeros(12,1);
xr = [0; 0; 0; 0; 0; 10; 0; 0; 0; 0; 0; 0];

% Prediction horizon
N = 10;

% Cast MPC problem to a QP: x = (x(0),x(1),...,x(N),u(0),...,u(N-1))
% - quadratic objective
P = blkdiag( kron(speye(N), Q), QN, kron(speye(N), R) );
% - linear objective
q = [repmat(-Q*xr, N, 1); -QN*xr; zeros(N*nu, 1)];
% - linear dynamics
Ax = kron(speye(N+1), -speye(nx)) + kron(sparse(diag(ones(N, 1), -1)), Ad);
Bu = kron([sparse(1, N); speye(N)], Bd);
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
nsim = 100;
buf=zeros(nsim,5);
for i = 1 : nsim
    % Solve
    res = prob.solve();

    % Check solver status
    if ~strcmp(res.info.status, 'solved')
        error('OSQP did not solve the problem!')
    end

    % Apply first control input to the plant
    ctrl = res.x((N+1)*nx+1:(N+1)*nx+nu);



    x0 = Ad*x0 + Bd*ctrl;
    buf(i,1)=x0(6);
    buf(i,2:5)=ctrl;

    % Update initial state
    l(1:nx) = -x0;
    u(1:nx) = -x0;
    prob.update('l', l, 'u', u);
    
    
end
plot(buf(:,1:5))
hold on




%%
syms k1 k2 c1 c2 J1 J2 u1_max u2_max u1_min u2_min q1_min q1_dot_min q2_min q2_dot_min tau1_min  tau2_min
syms q1_max q1_dot_max q2_max q2_dot_max tau1_max  tau2_max
syms q1_0 q1_dot_0 q2_0 q2_dot_0 tau1_0  tau2_0
syms q1_r q1_dot_r q2_r q2_dot_r tau1_r  tau2_r
syms Q1 Q2 Q3 Q4 Q5 Q6 R1 R2


N=1;
Ad = [  0       1       0       0    0    0;
    -k1/J1  -c1/J1   k2/J1     0   1/J1  0;
       0       0       0       1    0    0;
       0       0    -k2/J2  -c2/J2  0   1/J2;
       0       0       0       0    0    0;
       0       0       0       0    0    0];

Bd = [  0        0;
     k1/J1   -k2/J1;
       0        0;
       0      k2/J1;
       0        0;
       0        0];


[nx, nu] = size(Bd);

% Constraints
umin = [u1_min; u2_min];
umax = [u1_max; u2_max];

xmin = [q1_min; q1_dot_min; q2_min; q2_dot_min; tau1_min; tau2_min];
xmax = [q1_max; q1_dot_max; q2_max; q2_dot_max; tau1_max; tau2_max];

% Objective function
Q = diag([Q1 Q2 Q3 Q4 Q5 Q6]);
QN = Q;
R = diag([R1 R2]) ;

% Initial and reference states
x0 = [q1_0; q1_dot_0; q2_0; q2_dot_0; tau1_0; tau2_0];
xr = [q1_r; q1_dot_r; q2_r; q2_dot_r; tau1_r; tau2_r];

% Prediction horizon
% Cast MPC problem to a QP: x = (x(0),x(1),...,x(N),u(0),...,u(N-1))
% - quadratic objective
P = blkdiag( kron(speye(N), Q), QN, kron(speye(N), R) )
% - linear objective
q = [repmat(-Q*xr, N, 1); -QN*xr; zeros(N*nu, 1)]
% - linear dynamics
Ax = kron(speye(N+1), -speye(nx)) + kron(sparse(diag(ones(N, 1), -1)), Ad);
Bu = kron([sparse(1, N); speye(N)], Bd);
Aeq = [Ax, Bu];
leq = [-x0; zeros(N*nx, 1)];
ueq = leq;
% - input and state constraints
Aineq = speye((N+1)*nx + N*nu)
lineq = [repmat(xmin, N+1, 1); repmat(umin, N, 1)];
uineq = [repmat(xmax, N+1, 1); repmat(umax, N, 1)];
% - OSQP constraints
A = [Aeq; Aineq]
l = [leq; lineq]
u = [ueq; uineq]



%%
%model parameter setting 
J1 = 0.420*0.24215^2/3;
J2 = 0.225*0.1759^2/3;

k1 = 2*9.8;
k2 = 2*5.9;

c1 = 0.0001;
c2 = 0.0001;

A = [  0       1       0       0;
    -k1/J1  -c1/J1   k2/J1     0;
       0       0       0       1;
       0       0    -k2/J2  -c2/J2];

B = [  0       0;
     k1/J1  -k2/J1;
       0       0;
       0     k2/J1];

C = eye(size(A));

D = 0;

Wc = ctrb(A,B);
rankWc = rank(Wc)

% 제어 가능성 확인
if rankWc == size(A,1)
    disp('System is controllable.');
else
    disp('System is not controllable.');
end

Wo = obsv(A,C);
rankWo = rank(Wo)

% 관측 가능성 확인
if rankWo == size(A,1)
    disp('System is observable.');
else
    disp('System is not observable.');
end



