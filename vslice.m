N=10; % Number of timesteps


w = 
theta =
rho = 
Pi = 

% Time loop 
for n = 1:N
    % Initial guess at the our quantities, for the Newton solve
    % Could use a forward Euler step for the first guess.
    w = w_n;
    theta = theta_n;
    rho = rho_n;
    Pi = Pi_n;
    
    % Set the background reference states
    rho_ref = rho;
    theta_ref = theta;
    Pi_ref = Pi;
    
    % Calculate Rn_w, Rn_theta and Rn_rho. Sit on grid points.
    % Use a subfunction to calculate
    Rn_w = 
    Rn_theta = 
    Rn_rho = 
    
    % Set up for loop
    eta_d = eta_a - 
    % Departure point iteration
        for m = 1:M
        