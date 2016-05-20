% Script solving Burgers with a SL method, uniform grid.
clear all
clf

% Parameters
N=21;       % Level of spacial discretisation
tN = 200;
%Dt= 0.012;   % Timestep
theta_t=1/2;  % Theta for the Theta method in time.
theta_x=1/2;  % Theta for the Theta-method for departure points.
epsilon=0.002; % Epsilon in the PDE

x0 = 0;
x1 = 1;
Dx= (x1 - x0)./(N+1); % Space step

t0 = 0;
tmax = 1.5;
Dt = (tmax-t0)/tN;

% Initial and boundary conditions
u0 = @(x) sin(2*pi*x) + 1/2*sin(pi*x);
u_l = 0;
u_r = 0;

% Matrices to be used. Tridiagonal, so could be done more efficiently

% del2 - \delta_x^2, a tridiagonal matrix, finite difference second 
% derivative operator.
del2 = -2*eye(N) + diag(ones(N-1,1),1) + diag(ones(N-1,1),-1);
% LHS and RHS of SISL formulation, and the boundary term correction
M_RHS = eye(N) + Dt/(Dx^2) * (1 - theta_t) * epsilon * del2;
M_LHS = eye(N) - Dt/(Dx^2) * theta_t * epsilon * del2;
BC = [u_l/(Dx^2); zeros(N-2,1);u_r/(Dx^2)];

% Initialisation
X = x0 + (1:N)'*Dx;
Un = u0(X);

% For plotting...
tout = t0 + (0:tN)*Dt;
uout=zeros(length(tmax),N);
uout(1,:) = Un;
xout = X;
jj=1;

% Outer loop. Timestep
for ii = 1:(length(tout))
% Timestep initialisation
% Cubic spline of Un (including the constant end points)  
rhs0 = u_l + Dt/(Dx^2) * (1-theta_t) * epsilon *...
    (2*u_l - 5*Un(1) + 4*Un(2) - Un(3)); % This is an order h^2 approx'n to u''
rhsN = u_r + Dt/(Dx^2) * (1-theta_t) * epsilon *...
    (2*u_r - 5*Un(N) + 4*Un(N-1) - Un(N-2));
%pp_rhs = spline([0;X;1],[rhs0; (M_RHS * Un) ; rhsN]);
pp_rhs = interp1([x0;X;x1],...
  [rhs0; (M_RHS * Un + Dt*(1-theta_t)*epsilon*BC) ; rhsN],'linear','pp');
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
 U_A = M_LHS\(rhs_D + Dt*theta_t*epsilon*BC);
 % New guess at the departure points (theta-method) and recalculate 
 % RHS_D and U_A. 
 X_D = X - Dt*(theta_x*U_A + (1-theta_x)*Un);
 end
 rhs_D = ppval(pp_rhs, X_D);
 U_A = M_LHS\(rhs_D + Dt*theta_t*epsilon*BC);
 % End of Inner Loop
% Prep for next timestep.
Un = U_A;
plot([x0;X;x1],[u_l;U_A;u_r])
title(['t = ',num2str(tout(ii))]), ylim([-1 1.5]), drawnow()
jj=jj+1;
uout(jj,:)=Un;
end % for ii

% mk_video('test.avi','test',tout,uout,xout);
