clear
clc
Ts = 0.1;          % loop time [s]
d = 25.5*0.001;     %           [m]
sampling_num = 200; %Similation step
N = 10;              % Prediction horizon
    
J1 = 0.420*0.24215^2/3; %      [kgm2]
J2 = 0.225*0.1759^2/3;  %      [kgm2]
    
k1 = 15.05;             %      [N/m]
k2 = 7.396;             %      [N/m]
    
c1 = 0.01;              %      [Ns/m]
c2 = 0.01;              %      [Ns/m]
    
tau=0.0005;
delta_1=c1/tau-k1;
delta_2=c2/tau-k2;
x0 = [0; 0; 0; 0; 0; 0];                                 %initial state
x0_initial=x0;
    
q1=zeros(1,sampling_num);
q1_dot=zeros(1,sampling_num);
q2=zeros(1,sampling_num);
q2_dot=zeros(1,sampling_num);
phi1=zeros(1,sampling_num);
phi2=zeros(1,sampling_num);
    
theta1=zeros(1,sampling_num);
theta2=zeros(1,sampling_num);
    
X = [q1; q1_dot; q2; q2_dot; phi1; phi2];
X(:,1) = x0;
U = [theta1; theta2];

Xr = [q1; q1_dot; q2; q2_dot; phi1; phi2];
buf=zeros(1,N);

% Simulate in closed loop
for i = 1 : sampling_num
    %control value
    HZ=0.1; %hz
    w=HZ/(2*pi);
    xr = [0.5*sin(w*i); 0; sin(w*i); 0; 0; 0]; %rad rad/s rad rad/s rad rad
    xr = [deg2rad(20); 0; deg2rad(40); 0; 0; 0;]; %rad rad/s rad rad/s rad rad
    Xr(:,i+1) = xr;
    
    Ac = [      0           1          0         0            0            0;
          -(k1+k2)/J1  -(c1+c2)/J1   k2/J1     c2/J1          0            0;
                0           0          0         1            0            0;
              k2/J2       c2/J2     -k2/J2    -c2/J2          0            0;
                0           1          0         0            0            0;
                0           0          0         1            0            0];
    
    Bc = [   0       0   ;
           k1/J1  -k2/J1 ;
             0       0   ;
             0     k2/J2 ;
           -k1/J2     0   ;
             0     -k2/J2];

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
    xmin = [-(100.*3.141592/180.); -(100.*3.141592/180.); -(100.*3.141592/180.); -(100.*3.141592/180.); -(18.*3.141592/180.); -(18.*3.141592/180.)];
    xmax = [ (100.*3.141592/180.);  (100.*3.141592/180.);  (100.*3.141592/180.);  (100.*3.141592/180.);  (18.*3.141592/180.);  (18.*3.141592/180.)];
    
    % Objective function
    Q = diag([10 1 10 1 1 1]); %position control version
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
    
    
    for j=1:1:N
        delay_1=delta_1*exp(-Ts*j/tau);
        buf(1,j)=(k1+delay_1)/J1;
        delay_2=delta_2*exp(-Ts*j/tau);
    
        Bc_j = [         0                               0             ;
                  (k1+delay_1)/J1               -(k2+delay_2)/J1       ;
                         0                               0             ;
                         0                        (k2+delay_2)/J2      ;
                k1*exp(-Ts*j/tau)/(J2*tau)               0             ;
                         0                  k2*exp(-Ts*j/tau)/(J2*tau)];


                Bc_j = [         0                               0             ;
                  (k1+delay_1)/J1               -(k2+delay_2)/J1       ;
                         0                               0             ;
                         0                        (k2+delay_2)/J2      ;
                -k1/J2               0             ;
                         0                  -k2/J2];
    
        sys = ss(Ac, Bc_j, Cc, Dc);
        sysd = c2d(sys, Ts, 'zoh');
        Bd_j=sysd.B;
    
        A(j*nx+1:(j+1)*nx,nx*(N+1)+nu*(j-1)+1:nx*(N+1)+nu*j)=Bd_j;
    
    end
    
    % Create an OSQP object
    prob = osqp;
    
    % Setup workspace
    prob.setup(P, q, A, l, u);

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
X0 = x0_initial*ones(1, sampling_num); % Replicating X0 sampling_num times
t=linspace(0,sampling_num*Ts,sampling_num+1);
set(gcf, 'Color', 'w'); % Setting the background to white

% Drawing each subplot
subplot(4,2,1)
plot(t, rad2deg(X(1,:)), 'LineWidth', 2)
hold on
plot(0, rad2deg(X0(1,1)), 'bo', 'MarkerSize', 6, 'LineWidth', 1);
plot(t,rad2deg(Xr(1,:)), '--r', 'LineWidth', 3)
hold off
xlabel('t[s]', 'FontName', 'Times New Roman', 'FontSize', 12, 'FontWeight', 'bold')
ylabel('q1', 'FontName', 'Times New Roman', 'FontSize', 12, 'FontWeight', 'bold')

subplot(4,2,2)
plot(t,X(2,:), 'LineWidth', 2)
hold on
plot(0, rad2deg(X0(2,1)), 'bo', 'MarkerSize', 6, 'LineWidth', 1);
plot(t,rad2deg(Xr(2,:)), '--r', 'LineWidth', 2)
hold off
xlabel('t[s]', 'FontName', 'Times New Roman', 'FontSize', 12, 'FontWeight', 'bold')
ylabel('q1_dot', 'FontName', 'Times New Roman', 'FontSize', 12, 'FontWeight', 'bold')

subplot(4,2,3)
plot(t,rad2deg(X(3,:)), 'LineWidth', 2)
hold on
plot(0, rad2deg(X0(3,1)), 'bo', 'MarkerSize', 6, 'LineWidth', 1);
plot(t,rad2deg(Xr(3,:)), '--r', 'LineWidth', 2)
hold off
xlabel('t[s]', 'FontName', 'Times New Roman', 'FontSize', 12, 'FontWeight', 'bold')
ylabel('q2', 'FontName', 'Times New Roman', 'FontSize', 12, 'FontWeight', 'bold')

subplot(4,2,4)
plot(t,rad2deg(X(4,:)), 'LineWidth', 2)
hold on
plot(0, rad2deg(X0(4,1)), 'bo', 'MarkerSize', 6, 'LineWidth', 1);
plot(t,rad2deg(Xr(4,:)), '--r', 'LineWidth', 2)
hold off
xlabel('t[s]', 'FontName', 'Times New Roman', 'FontSize', 12, 'FontWeight', 'bold')
ylabel('q2_dot', 'FontName', 'Times New Roman', 'FontSize', 12, 'FontWeight', 'bold')

subplot(4,2,5)
plot(t,rad2deg(X(5,:)), 'LineWidth', 2)
hold on
plot(0, rad2deg(X0(5,1)), 'bo', 'MarkerSize', 6, 'LineWidth', 1);
plot(t,rad2deg(Xr(5,:)), '--r', 'LineWidth', 2)
hold off
xlabel('t[s]', 'FontName', 'Times New Roman', 'FontSize', 12, 'FontWeight', 'bold')
ylabel('phi1', 'FontName', 'Times New Roman', 'FontSize', 12, 'FontWeight', 'bold')

subplot(4,2,6)
plot(t,rad2deg(X(6,:)), 'LineWidth', 2)
hold on
plot(0, rad2deg(X0(6,1)), 'bo', 'MarkerSize', 6, 'LineWidth', 1);
plot(t,rad2deg(Xr(6,:)), '--r', 'LineWidth', 2)
hold off
xlabel('t[s]', 'FontName', 'Times New Roman', 'FontSize', 12, 'FontWeight', 'bold')
ylabel('phi2', 'FontName', 'Times New Roman', 'FontSize', 12, 'FontWeight', 'bold')

subplot(4,2,7)
plot(t,rad2deg(U(1,:)), 'LineWidth', 2)
xlabel('t[s]', 'FontName', 'Times New Roman', 'FontSize', 12, 'FontWeight', 'bold')
ylabel('theta1', 'FontName', 'Times New Roman', 'FontSize', 12, 'FontWeight', 'bold')

subplot(4,2,8)
plot(t,rad2deg(U(2,:)), 'LineWidth', 2)
xlabel('t[s]', 'FontName', 'Times New Roman', 'FontSize', 12, 'FontWeight', 'bold')
ylabel('theta2', 'FontName', 'Times New Roman', 'FontSize', 12, 'FontWeight', 'bold')


%%
%model parameter setting 
J1 = 0.420*0.24215^2/3;
J2 = 0.225*0.1759^2/3;

k1 = 2*9.8;
k2 = 2*5.9;

c1 = 0.0001;
c2 = 0.0001;

    A = [      0           1          0         0            0            0;
          -(k1+k2)/J1  -(c1+c2)/J1   k2/J1     c2/J1          0            0;
                0           0          0         1            0            0;
              k2/J2       c2/J2     -k2/J2    -c2/J2          0            0;
                0           1          0         0            0            0;
                0           0          0         1            0            0];

    B = [   0       0   ;
           k1/J1  -k2/J1 ;
             0       0   ;
             0     k2/J2 ;
           k1/J2     0   ;
             0     k2/J2];

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
