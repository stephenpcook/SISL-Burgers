% Solving Burgers with a SL method.
clear all
clf

% Parameters
N=21;       % Level of spacial discretisation
Dx= 1./(N+1); % Space step
tN = 200;
%Dt= 0.012;   % Timestep
theta_t=1/2;  % Theta for the Theta method in time.
theta_x=1/2;  % Theta for the Theta-method for departure points.
epsilon=0.002; % Epsilon in the PDE
tmax = 1.5;

Dt = tmax/tN;


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

% Initialisation
X = (1:N)'./(N+1);
Un = u0(X);

% For plotting...
tout=(0:Dt:tmax)';
uout=zeros(length(tmax),N);
uout(1,:) = Un;
xout = X;
jj=1;

% Outer loop. Timestep
for t = 0:Dt:(tmax-Dt)
% Timestep initialisation
% Cubic spline of Un (including the constant end points)  
rhs0 = -Dt/(Dx^2) * (1-theta_t) * epsilon *...
    (5*Un(1) + 4*Un(2) - Un(3)); % This is an order h^2 approx'n to u''
rhsN = -Dt/(Dx^2) * (1-theta_t) * epsilon *...
    (5*Un(N) + 4*Un(N-1) - Un(N-2));
%pp_rhs = spline([0;X;1],[rhs0; (M_RHS * Un) ; rhsN]);
pp_rhs = interp1([0;X;1],[rhs0; (M_RHS * Un) ; rhsN],'linear','pp');
%   TODO
%   This is actually wrong at the moment. Want to take into account
%   the second derivative term at the left and right points.
%   Since this is in the known LHS, could make this higher order;
%   doesn't have to lead to a tridiagonal matrix to solve with.

% Initial guess of departure points.
X_D = X - Dt*Un;
 % Inner loop. 
 for k = 1:2
 % Evaluate RHS at the departure points.
 rhs_D = ppval(pp_rhs, X_D);
 % Solve the implicit equation (Thomas algorithm implemented later)
 U_A = M_LHS\rhs_D;
 % New guess at the departure points (theta-method) and recalculate 
 % RHS_D and U_A. 
 X_D = X - Dt*(theta_x*U_A + (1-theta_x)*Un);
 rhs_D = ppval(pp_rhs, X_D);
 end
 U_A = M_LHS\rhs_D;
 % End of Inner Loop
% Prep for next timestep.
Un = U_A;
%plot([u_l;U_A;u_r]), title(['t = ',num2str(t)]), ylim([-1 1.5]), drawnow()
jj=jj+1;
uout(jj,:)=Un;
end % for t

% mk_video('test.avi','test',tout,uout,xout);


