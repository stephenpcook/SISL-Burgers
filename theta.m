% Solving Burgers with theta method.
clear all
clf

% Parameters
N= 50;       % Level of spacial discretisation
Dx= 1./(N+1); % Space step
Dt= 0.012;   % Timestep
theta_t=0;  % Theta for the Theta method in time.
theta_x=1/2;  % Theta for the Theta-method for departure points.
epsilon=0; % Epsilon in the PDE
tmax = 10;


% Initial and boundary conditions
u0 = @(x) sin(2*pi*x) + 1/2*sin(pi*x);
u_l = 0;
u_r = 0;

% Matrices to be used. Tridiagonal, so could be done more efficiently

% del2 - \delta_x^2, a tridiagonal matrix, finite difference second 
% derivative operator.
del2 = 2*eye(N) - diag(ones(N-1,1),1) - diag(ones(N-1,1),-1);
% LHS and RHS of SISL formulation
M_RHS = eye(N) - Dt/(Dx^2) * (1 - theta_t) * epsilon * del2;
M_LHS = eye(N) + Dt/(Dx^2) * theta_t * epsilon * del2;

F_LHS = @(u) Dt/Dx*theta_t/2*[(2)^2-u(1)^2;diff(u.^2)];
F_RHS = @(u) - Dt/Dx*(1-theta_t)/2*[u(2)^2-u(1)^2;diff(u.^2)];

% Initialisation
X = (1:N)'./(N+1);
Un = u0(X);

% Outer loop. Timestep
for t = 0:Dt:tmax
% Timestep initialisation

u_rhs = (M_RHS * Un + F_RHS(Un) - F_LHS(Un));

U_A = M_LHS\u_rhs;
 % End of Inner Loop
% Prep for next timestep.
Un = U_A;
plot([u_l;U_A;u_r]), title(['t = ',num2str(t)]), ylim([-1 1.5]), drawnow()
end % for t