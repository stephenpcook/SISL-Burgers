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
x0 = -1;
x1 = 4;
%Dx= (x1-x0)./(N+1); % Space step

t0 = 0;
tmax = 1.5;   % Final time
%Dt = (tmax-t0)/tN;

% Initial and boundary conditions
%u0 = @(x) (sin(2*pi*x) + 1/2*sin(pi*x)) + 1;
%u_l = 1; u_r = 1;
c = 1;
alpha_0 = 0.1;
u0 = @(x) c - alpha_0*tanh(alpha_0/(2*epsilon)*(x - c*t0));
u_l = u0(x0);
u_r = u0(x1);

% Plot parameters
%plotlims = [-0.5, 1.5];
plotlims = [c-alpha_0-0.1, c+alpha_0+0.1];
plotting = 1;

% Mesh Parameters
%mesh = 'static';
%mesh = 'prescribed';
mesh_movement = 'moving-exact';  % Mesh movement type
%mesh = 'moving-relax';  % Mesh movement type
limiter = 1;        % Flux limiter for interpolation
interpolation = 'linear';
%interpolation = 'CSpline';
%interpolation = 'CLagrange';
%interpolation = 'ENO';
%interpolation = 'pchip';

% Monitor function parameters
b = 0.1;
m=@(x,u,uprime) sqrt(b + uprime.^2);
%m=@(x,u,uprime)ones(size(u))
p_smooth = 5;
tau =1;
with_euler = 1;

save('params_default',...
  ...% General Parameters
  'K','theta_t','theta_x','epsilon',...
  ...% Domain Parameters
  'x0', 'x1', 't0', 'tmax',...
  ...% Initial and boundary conditions
  'c', 'alpha_0', 'u0', 'u_l', 'u_r', ...
  ...% Plot parameters
  'plotlims', 'plotting', ...
  ...% Mesh Parameters
  'mesh_movement', 'b', 'm', 'p_smooth', 'tau', 'with_euler',...
  ...% Interpolation Parameters
  'limiter', 'interpolation'...
  );
end % function gen_defaults
