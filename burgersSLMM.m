% Solving Burgers with a SL method.
clear all
clf

% Parameters
N=20;       % Level of spacial discretisation
Dt= 0.012;   % Timestep
tmax = 2;    % Final time
theta_t=1/2;  % Theta for the Theta method in time.
theta_x=1/2;  % Theta for the Theta-method for departure points.
epsilon=0.01; % Epsilon in the PDE
% Mesh Parameters
mesh = 'static';
%mesh = 'perscribed';
mesh = 'moving';
b = 0.1;

Dx= 1./(N+1); % Space step


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
%X = X.^2;
Un = u0(X);
X_An = X;
TT = 0:Dt:tmax;
XX = zeros(length(TT),length(X));
XX(1,:) = X;

del2n1 = eye(N);
del2n = eye(N);

% Outer loop. Timestep
for tt = 1:length(TT)
    t = TT(tt);
% Timestep initialisation
% Cubic spline of Un (including the constant end points)  

% Mesh movement
switch mesh
    case 'static'
        X_An1 = X;
    case 'perscribed'
        % Perscribed mesh, loaded from file innit.
        X_An1 = X;
    case 'moving'
        monitor.type = 'points';
        monitor.x = [0;X_An;1];
        %monitor.M = sqrt(0.1 + [u_l;Un;u_r].^2);
        DUDX = [u_l;diff([u_l;Un;u_r])]./[1;diff([0;X_An;1])];
        monitor.M = 1/2+1/2*sqrt(b + DUDX.^2);
        X_An1 = Eqd1dExact([0;X_An;1], monitor);
        X_An1 = X_An1(2:end-1);
 %       X_An1 = (X_An1 + X_An)./2;
    otherwise
        disp('Unknown movement type, using static')
        X_An1 = X;
end % switch

DX_An = diff([0;X_An;1]);
DX_An1 = diff([0;X_An1;1]);

del2n = del2n1;

% Build up the LHS and RHS tridiagonal matrix. Can recycle these 
% and get 2 uses out of each del2 matrix (only redo one)
del2n1(1,1) = -2/(DX_An1(2)*DX_An1(1));
del2n1(1,2) = 2/(DX_An1(2)*(DX_An1(2)+DX_An1(1)));
for ii=2:(N-1)
    del2n1(ii,ii-1) = 2/(DX_An1(ii)*(DX_An1(ii+1)+DX_An1(ii)));
    del2n1(ii,ii) = -2/(DX_An1(ii)*DX_An1(ii+1));
    del2n1(ii,ii+1) = 2/(DX_An1(ii+1)*(DX_An1(ii+1)+DX_An1(ii)));
end
del2n1(N,N-1) = 2/(DX_An1(N-1)*(DX_An1(N)+DX_An1(N-1)));
del2n1(N,N) = -2/(DX_An1(N)*DX_An1(N-1));
M_LHS = eye(N) - Dt * theta_t * epsilon * del2n1;

%del2n(1,1) = -2/(DX_An(2)*DX_An(1));
%del2n(1,2) = 2/(DX_An(2)*(DX_An(2)+DX_An(1)));
%for ii=2:(N-1)
%    del2n(ii,ii-1) = 2/(DX_An(ii)*(DX_An(ii+1)+DX_An(ii)));
%    del2n(ii,ii) = -2/(DX_An(ii)*DX_An(ii+1));
%    del2n(ii,ii+1) = 2/(DX_An(ii+1)*(DX_An(ii+1)+DX_An(ii)));
%end
%del2n(N,N-1) = 2/(DX_An(N-1)*(DX_An(N)+DX_An(N-1)));
%del2n(N,N) = -2/(DX_An(N)*DX_An(N-1));
M_RHS = eye(N) + Dt * (1 - theta_t) * epsilon * del2n;

% Make the spline for interpolation
rhs0 = -Dt/(Dx^2) * (1-theta_t) * epsilon *...
    (5*Un(1) + 4*Un(2) - Un(3)); % This is an order h^2 approx'n to u''
rhsN = -Dt/(Dx^2) * (1-theta_t) * epsilon *...
    (5*Un(N) + 4*Un(N-1) - Un(N-2));
pp_rhs = spline([0;X_An;1],[rhs0; (M_RHS * Un) ; rhsN]);

% Initial guess of departure points.
X_D = X_An1 - Dt*Un;
 % Inner loop. 
 for k = 1:2
 % Evaluate RHS at the departure points.
 rhs_D = ppval(pp_rhs, X_D);
 % Solve the implicit equation (Thomas algorithm implemented later)
 U_A = M_LHS\rhs_D;
 % New guess at the departure points (theta-method) and recalculate 
 % RHS_D and U_A. 
 X_D = X_An1 - Dt*(theta_x*U_A + (1-theta_x)*Un);
 rhs_D = ppval(pp_rhs, X_D);
 end
 U_A = M_LHS\rhs_D;
 % End of Inner Loop


%%% Prep for next timestep. %%%
X_An = X_An1;

Un = U_A;
%%% And plot %%% 
plot([0;X_An1;1],[u_l;U_A;u_r]), title(['t = ',num2str(t)]), ylim([-1 1.5])
hold on
plot(monitor.x,monitor.M,'g-')
hold off
drawnow()
XX(tt,:) = X_An;
bigX_D(tt,:) = X_D;
end % for t
pause
for i = 1:N
    plot(XX(:,i),TT)
    hold on
    plot(bigX_D(:,i),TT,'g-')
end