clear
clc
%envirment setting
Ts = 0.01; %s
d = 25.5*0.001;
sampling_num = 5000;
N = 10; % Prediction horizon

J1 = 0.420*0.24215^2/3;
J2 = 0.225*0.1759^2/3;

k1 = 2*9.8;
k2 = 2*5.9;

c1 = 0.1;
c2 = 0.1;

%control value
xr = [deg2rad(20); 0; deg2rad(20); 0; deg2rad(3); 0]; %rad rad/2 rad rad/2 Nm Nm
x0 = [0; 0; 0; 0; 0; 0]; %initial state
x0_initial=x0;

q1=zeros(1,sampling_num);
q1_dot=zeros(1,sampling_num);
q2=zeros(1,sampling_num);
q2_dot=zeros(1,sampling_num);
phi1=zeros(1,sampling_num);
phi2=zeros(1,sampling_num);

theta1=zeros(1,sampling_num);
theta1_dot=zeros(1,sampling_num);
theta2=zeros(1,sampling_num);
theta2_dot=zeros(1,sampling_num);

theta1_dot_Lim=zeros(1,sampling_num);
theta2_dot_Lim=zeros(1,sampling_num);

X = [q1; q1_dot; q2; q2_dot; phi1; phi2];
X(:,1) = x0;
U = [theta1; theta1_dot; theta2; theta2_dot];
U_Lim = [theta1_dot_Lim; theta2_dot_Lim];

A = [  0       1       0       0       0       0;
    -k1/J1  -c1/J1   k2/J1   c2/J1   k1/J1     0;
       0       0       0       1       0       0;
       0       0    -k2/J2  -c2/J2     0     k2/J2;
       0       1       0       0       0       0;
       0       0       0       1       0       0];

B = [  0       1       0       0;
     k1/J1   c1/J1  -k2/J1  -c2/J1;
       0       0       0       1;
       0       0     k2/J2   c2/J2;
       0      -1       0       0;
       0       0       0      -1];

C = eye(size(A));

D = 0;
sys = ss(A, B, C, D);
sysd = c2d(sys, Ts, 'zoh');
Ad=sysd.A;
Bd=sysd.B;
[nx, nu] = size(Bd);

% Constraints_py ver
umin = [-deg2rad(100)-d*deg2rad(100); 0; -deg2rad(100)-d*deg2rad(100); 0];
umax = [ deg2rad(100)+d*deg2rad(100); 0;  deg2rad(100)-d*deg2rad(100); 0];
xmin = [-deg2rad(100); -deg2rad(100); -deg2rad(100); -deg2rad(100); -deg2rad(10); -deg2rad(10)];
xmax = [ deg2rad(100);  deg2rad(100);  deg2rad(100);  deg2rad(100);  deg2rad(10);  deg2rad(10)];

% Objective function
Q = diag([10 1 10 1 30 10]); %position control version
QN = Q;
R = diag([0.05 0 0.05 0]);

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
for i = 1 : sampling_num
    % Solve
    res = prob.solve();
    
    
    % Apply first control input to the plant
    ctrl = res.x((N+1)*nx+1:(N+1)*nx+nu);
    
    

    x0 = Ad*x0 + Bd*(ctrl);

    X(:,i+1) = x0;
    U(:,i+1) = ctrl;
    

    % Update initial state
    l(1:nx) = -x0;
    u(1:nx) = -x0;
    theta1_dot_limit = (U(2,i+1)-U(2,i))/2;
    theta2_dot_limit = (U(4,i+1)-U(4,i))/2;
    U_Lim(:,i+1)=[theta1_dot_limit, theta1_dot_limit];


    l_new=[l(2*nx*(N+1)+1); theta1_dot_limit;  l(2*nx*(N+1)+1+2);  theta2_dot_limit];
    l(2*nx*(N+1)+1:end)=repmat(l_new,N,1);

    u_new=[u(2*nx*(N+1)+1); theta1_dot_limit;  u(2*nx*(N+1)+1+2);  theta2_dot_limit];
    u(2*nx*(N+1)+1:end)=repmat(u_new,N,1);

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
plot(Xr(5,:), '--r', 'LineWidth', 2)
hold off
ylim([-1 1])
xlabel('step', 'FontName', 'Times New Roman', 'FontSize', 12, 'FontWeight', 'bold')
ylabel('phi1', 'FontName', 'Times New Roman', 'FontSize', 12, 'FontWeight', 'bold')

subplot(5,2,6)
plot(X(6,:), 'LineWidth', 2)
hold on
plot(1, X0(6,1), 'bo', 'MarkerSize', 6, 'LineWidth', 1);
plot(Xr(6,:), '--r', 'LineWidth', 2)
hold off
ylim([-1 1])
xlabel('step', 'FontName', 'Times New Roman', 'FontSize', 12, 'FontWeight', 'bold')
ylabel('phi2', 'FontName', 'Times New Roman', 'FontSize', 12, 'FontWeight', 'bold')

subplot(5,2,7)
plot(U(1,:), 'LineWidth', 2)
ylim([-1 1])
xlabel('step', 'FontName', 'Times New Roman', 'FontSize', 12, 'FontWeight', 'bold')
ylabel('theta1', 'FontName', 'Times New Roman', 'FontSize', 12, 'FontWeight', 'bold')

subplot(5,2,8)
plot(U(3,:), 'LineWidth', 2)
ylim([-1 1])
xlabel('step', 'FontName', 'Times New Roman', 'FontSize', 12, 'FontWeight', 'bold')
ylabel('theta2', 'FontName', 'Times New Roman', 'FontSize', 12, 'FontWeight', 'bold')

subplot(5,2,9)
plot(U(2,:), 'LineWidth', 2)
ylim([-1 1])
hold on
plot(U_Lim(1,:), '--r', 'LineWidth', 2)
xlabel('step', 'FontName', 'Times New Roman', 'FontSize', 12, 'FontWeight', 'bold')
ylabel('theta1_dot', 'FontName', 'Times New Roman', 'FontSize', 12, 'FontWeight', 'bold')

subplot(5,2,10)
plot(U(4,:), 'LineWidth', 2)
ylim([-1 1])
hold on
plot(U_Lim(2,:), '--r', 'LineWidth', 2)
xlabel('step', 'FontName', 'Times New Roman', 'FontSize', 12, 'FontWeight', 'bold')
ylabel('theta2_dot', 'FontName', 'Times New Roman', 'FontSize', 12, 'FontWeight', 'bold')




%%
%model parameter setting 
J1 = 0.420*0.24215^2/3;
J2 = 0.225*0.1759^2/3;

k1 = 2*9.8;
k2 = 2*5.9;

c1 = 0.0001;
c2 = 0.0001;

A = [  0       1       0       0       0       0;
    -k1/J1  -c1/J1   k2/J1   c2/J1   k1/J1     0;
       0       0       0       1       0       0;
       0       0    -k2/J2  -c2/J2     0     k2/J2;
       0       1       0       0       0       0;
       0       0       0       1       0       0];

B = [  0       1       0       0;
     k1/J1   c1/J1  -k2/J1  -c2/J1;
       0       0       0       1;
       0       0     k2/J2   c2/J2;
       0      -1       0       0;
       0       0       0      -1];

C = eye(size(A));

D = 0;


Wc = ctrb(A,B);
rankWc = rank(Wc);

[U, S, V] = svd(Wc);
singular_values = diag(S)

% 제어 가능성 확인
% if rankWc == size(A,1)
%     disp('System is controllable.');
% else
%     disp('System is not controllable.');
% end
% 
% Wo = obsv(A,C);
% rankWo = rank(Wo)
% 
% % 관측 가능성 확인
% if rankWo == size(A,1)
%     disp('System is observable.');
% else
%     disp('System is not observable.');
% end



