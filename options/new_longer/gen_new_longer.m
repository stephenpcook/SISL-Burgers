function gen_new_longer()
% Creates param_defaults.mat

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%               Parameters              %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% General Parameters
K=10;         % Departure point iterations
epsilon=1e-4; % Epsilon in the PDE

% Domain Parameters
x_l = -5;
x_r = 5;

t0 = 0;
tmax = 2;   % Final time

% Initial and boundary conditions
c = 1;
alpha_0 = 0.1;
u0 = @(x) c - alpha_0*tanh(alpha_0/(2*epsilon)*(x - c*t0));
u_l = u0(x_l);
u_r = u0(x_r);

% Plot parameters
plotting = 0;
plotlims = [c-alpha_0-0.1, c+alpha_0+0.1];

%% new_longer_static

% Mesh Parameters
mesh_movement = 'static';
interpolation = 'linear';

save('new_longer_static'...
  ...% General Parameters
  ,'K', 'epsilon'...
  ...% Domain Parameters
  ,'x_l', 'x_r', 't0', 'tmax'...
  ...% Initial and boundary conditions
  ,'c', 'alpha_0', 'u0', 'u_l', 'u_r'...
  ...% Plot parameters
  ,'plotting', 'plotlims'...
  );

%% new_longer_static_hermite

% Mesh Parameters
mesh_movement = 'static';
interpolation = 'hermite';
limiter = 0;

save('new_longer_static_hermite'...
  ...% General Parameters
  ,'K', 'epsilon'...
  ...% Domain Parameters
  ,'x_l', 'x_r', 't0', 'tmax'...
  ...% Initial and boundary conditions
  ,'c', 'alpha_0', 'u0', 'u_l', 'u_r'...
  ...% Plot parameters
  ,'plotting', 'plotlims'...
  );


end % function gen_new_longer
