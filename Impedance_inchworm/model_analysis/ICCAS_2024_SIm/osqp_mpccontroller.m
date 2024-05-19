function ctrl = osqp_mpccontroller(xr_and_x0)

xr=xr_and_x0(1:8);
x0=xr_and_x0(9:end);


% MPC parameter
umin = [-deg2rad(100) ; -deg2rad(100)];
umax = [ deg2rad(100) ;  deg2rad(100)];
xmin = [-deg2rad(100) ; -deg2rad(100) ; -deg2rad(100) ; -deg2rad(100) ; umin ; umin];
xmax = [ deg2rad(100) ;  deg2rad(100) ;  deg2rad(100) ;  deg2rad(100) ; umax ; umax];

% envirment setting
Ts = 0.01; % loop rate [s]
N = 10;    % Prediction horizon

% Bessel Filter state space
w0_1=30; % cut-off fre [rad/s]
w0_2=30; % cut-off fre [rad/s]

% SEA state space
J1 = 0.420*0.24215^2/12+0.420*0.13475^2;
J2 = 0.225*0.1759^2/12+0.225*0.0562^2;
k1 = 15.05;
k2 = 7.396;
c1 = 0.1;
c2 = 0.1;

A_SEA = [     0               1           0          0;
         -(k1+k2)/J1    -(c1+c2)/J1     k2/J1      c2/J1;
              0               0           0          1;
            k2/J2           c2/J2      -k2/J2     -c2/J2];

B_SEA = [  0      0        0         0;
         k1/J1  c1/J1   -k2/J1    -c2/J1;
           0      0        0         0;
           0      0      k2/J2    c2/J2];


A_motor=[     0        1        0           0;
         -3*w0_1^2  -3*w0_1     0           0;
              0        0        0           1;
              0        0    -3*w0_2^2    -3*w0_2];

B_motor=[    0       0;
          3*w0_1^2   0;
             0       0;
             0    3*w0_2^2];


% Objective function
Q = diag([100 1 10 1 0 0 0 0]);       %[q1 q1_dot q2 q2_dot theta1 theta1_dot theta2 theta2_dot]
QN = Q;                              %[q1 q1_dot q2 q2_dot theta1 theta1_dot theta2 theta2_dot]

R_theta_des = diag([0.01 0.01]);         %[theta1_des theta2_des]
R = R_theta_des;                     %[theta1_des theta2_des]

A_model=[       A_SEA                  B_SEA;
                zeros(size(A_motor))   A_motor];

B_model=[       zeros(size(B_SEA,1),2); 
                B_motor              ];

[nx, nu] = size(B_model);
C_model=eye(nx);

D_model=zeros(size(B_model));

sys = ss(A_model, B_model, C_model, D_model);
sysd = c2d(sys, Ts, 'zoh');
Ad_model=sysd.A;
Bd_model=sysd.B;

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
res = prob.solve();
ctrl = res.x((N+1)*nx+1:(N+1)*nx+nu);



end
