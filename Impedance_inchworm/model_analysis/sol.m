clear
clc
%envirment setting
Ts = 0.1;          % loop time [s]
d = 25.5*0.001;     %           [m]
sampling_num = 20; %Similation step
N = 5;              % Prediction horizon

J1 = 0.420*0.24215^2/3; %      [kgm2]
J2 = 0.225*0.1759^2/3;  %      [kgm2]

k1 = 2*9.8;             %      [N/m]
k2 = 2*5.9;             %      [N/m]

c1 = 0.01;              %      [Ns/m]
c2 = 0.01;              %      [Ns/m]

tau=0.0005;
delta_1=c1/tau-k1;
delta_2=c2/tau-k2;


%control value
xr = [deg2rad(30); 0; deg2rad(20)+deg2rad(20); 0]; %rad rad/s rad rad/s
x0 = [0.3; 0; 0.3; 0];                                 %initial state
x0_initial=x0;

q1=zeros(1,sampling_num);
q1_dot=zeros(1,sampling_num);
q2=zeros(1,sampling_num);
q2_dot=zeros(1,sampling_num);

theta1=zeros(1,sampling_num);
theta2=zeros(1,sampling_num);

X = [q1; q1_dot; q2; q2_dot];
X(:,1) = x0;
U = [theta1; theta2];

Ac = [      0           1          0         0   ;
      -(k1+k2)/J1  -(c1+c2)/J1   k2/J1     c2/J1 ;
            0           0          0         1   ;
          k2/J2       c2/J2     -k2/J2    -c2/J2];

Bc = [   0       0   ;
       k1/J1  -k2/J1 ;
         0       0   ;
         0     k2/J2];

Cc = eye(size(Ac));

Dc = 0;
sys = ss(Ac, Bc, Cc, Dc);
sysd = c2d(sys, Ts, 'zoh');
Ad=sysd.A;
Bd=sysd.B;
[nx, nu] = size(Bd);

% Constraints_py ver
umin = [-(100.*3.141592/180.)-d*(100.*3.141592/180.); -(100.*3.141592/180.)-d*(100.*3.141592/180.)];
umax = [ (100.*3.141592/180.)+d*(100.*3.141592/180.);  (100.*3.141592/180.)+d*(100.*3.141592/180.)];
xmin = [-(100.*3.141592/180.); -(90.*3.141592/180.); -(100.*3.141592/180.); -(90.*3.141592/180.);];
xmax = [(100.*3.141592/180.); (90.*3.141592/180.); (100.*3.141592/180.); (90.*3.141592/180.);];

% Objective function
Q = diag([10 1 10 1]); %position control version
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

buf=zeros(1,N);
for i=1:1:N
    delay_1=delta_1*exp(-Ts*i/tau);
    buf(1,i)=(k1+delay_1)/J1;
    delay_2=delta_2*exp(-Ts*i/tau);

    Bc_i = [      0                  0       ;
           (k1+delay_1)/J1  -(k2+delay_2)/J1 ;
                  0                  0       ;
                  0          (k2+delay_2)/J2];

    sys = ss(Ac, Bc_i, Cc, Dc);
    sysd = c2d(sys, Ts, 'zoh');
    Bd_i=sysd.B;

    A(i*nx+1:(i+1)*nx,nx*(N+1)+nu*(i-1)+1:nx*(N+1)+nu*i)=Bd_i;

end
plot(buf)

% Create an OSQP object
prob = osqp;

% Setup workspace
prob.setup(P, q, A, l, u);


%
% Simulate in closed loop
for i = 1 : sampling_num
    dw=i
    res = prob.solve();
    if ~strcmp(res.info.status, 'solved')
        error('OSQP did not solve the problem!')
    end
    ctrl = res.x((N+1)*nx+1:(N+1)*nx+nu);

    x0 = Ad*x0 + Bd*(ctrl);
    % x0 = Ad*x0;
    X(:,i+1) = x0;
    U(:,i+1) = ctrl;


    % Update initial state
    l(1:nx) = -x0;
    u(1:nx) = -x0;
    prob.update('l', l, 'u', u);


end

% Assuming data preparation for Xr and X0
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
