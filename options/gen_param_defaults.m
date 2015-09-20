function gen_param_defaults()
% Creates param_defaults.mat

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%               Parameters              %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% General Parameters
K=10;         % Departure point iterations
theta_t=1/2;  % Theta for the Theta method in time.
theta_x=1/2;  % Theta for the Theta-method for departure points.
epsilon=0.0001 ; % Epsilon in the PDE

% Domain Parameters
x_l = -1;
x_r = 4;
%Dx= (x_r-x_l)./(N+1); % Space step

t0 = 0;
tmax = 1.5;   % Final time
%Dt = (tmax-t0)/tN;

% Initial and boundary conditions
%u0 = @(x) (sin(2*pi*x) + 1/2*sin(pi*x)) + 1;
%u_l = 1; u_r = 1;
c = 1;
alpha_0 = 0.1;
u0 = @(x) c - alpha_0*tanh(alpha_0/(2*epsilon)*(x - c*t0));
u_l = u0(x_l);
u_r = u0(x_r);

% Plot parameters
%plotlims = [-0.5, 1.5];
plotting = 1;
plotlims = [c-alpha_0-0.1, c+alpha_0+0.1];

% Mesh Parameters
mesh_movement = 'static';
%mesh_movement = 'prescribed';
%mesh_movement = 'moving-exact';  % Mesh movement type
%mesh_movement = 'moving-relax';  % Mesh movement type
interpolation = 'linear';
%interpolation = 'CSpline';
%interpolation = 'CLagrange';
%interpolation = 'ENO';
%interpolation = 'pchip';
limiter = 1;        % Flux limiter for interpolation

% Monitor function parameters
b = 0.1;
m=@(x,u,uprime) sqrt(b + uprime.^2);
%m=@(x,u,uprime)ones(size(u))
p_smooth = 5;
tau =1;
with_euler = 1;

save('params_default'...
  ...% General Parameters
  ,'K','theta_t','theta_x','epsilon'...
  ...% Domain Parameters
  ,'x_l', 'x_r', 't0', 'tmax'...
  ...% Initial and boundary conditions
  ,'c', 'alpha_0', 'u0', 'u_l', 'u_r'...
  ...% Plot parameters
  ,'plotting', 'plotlims'...
  ...% Mesh Parameters
  ,'mesh_movement'...
  ...% Monitor function parameters
  ...%,'m', 'p_smooth', 'tau', 'with_euler'...
  ...% Interpolation Parameters
  ,'interpolation', 'limiter'...
  );

end % function gen_defaults
