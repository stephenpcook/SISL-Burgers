function [Un,X_A,bigXstar] = burg2(N, Nt, param_file)
%BURG2 A light-weight SISL burgers solver with uniform grids and interp1
%
% [Un,X_A] = BURG2(N,Nt) returns the end-time numerical solution and its
% associated uniform mesh with parameters defined by the options file
% options/params_default.mat, defined by gen_param_defaults.
%
% [Un,X_A] = BURG2(N,Nt,param_file) is as above, but takes parameters from
% the .mat parameter files from param_file. Only takes domain and plotting
% options, and ignores any options for interpolation or mesh movement.
%
% [Un,X_A,bigXstar] = BURG2(N,Nt, ... ) tracks the front, making a call to
% get_m_x finding the points x_star(t) such that
%
%   U(x_star) = (u(x_l) + u(x_r))/2.
%
% See also: BURGERSSLMM GEN_PARAM_DEFAULTS GET_M_X

%% General Parameters
%%N = 300;
%%Nt = 3000;
%K = 40;
%theta_x = 0.5;
%epsilon = 0.0001;
%
%c = 1;
%alpha_0 = 0.1;
%
%% Domain parameters
%x_l = -1;
%x_r = 4;
%
%t0 = 0;
%tmax = 1.5;
%
%% Initial and Boundary Conditions
%
%%x_l = 0;
%%x_r = 1;
%%u0 = @(x) sin(2*pi*x) + 1/2*sin(pi*x) - 0.5;
%
%u0 = @(x) c - alpha_0*tanh(alpha_0/(2*epsilon)*(x - c*t0));
%u_l = u0(x_l);
%u_r = u0(x_r);
%
%%tmax = 15;
%%x_l = -10;
%%x_r = 40;
%
%% Plot parameters
%%plotlims = [-2.5 2.5];
%plotting = 0;
%plotlims = [0 1.6];

old_path = addpath([pwd,'\options']);

load('params_default',...
  'K', 'theta_t', 'theta_x', 'epsilon', ...
  'x_l', 'x_r', 't0', 'tmax',...
  'c', 'alpha_0', 'u0', 'u_l', 'u_r',...
  'plotting', 'plotlims');
if nargin==3
  load(param_file);
end

if nargout==3
  track_front = 1;
else
  track_front = 0;
end

tmax = 3;
plotting = 1;
plottimes = [floor(40/3), ceil(40*2/3), 40, floor(40*4/3)];

% Setup
Dx = (x_r-x_l)/(N+1);
Dt = (tmax-t0)/Nt;

X_A_bc = (x_l:Dx:x_r)';
X_A = X_A_bc(2:end-1);

delta2 = 1/(Dx^2) * ...
    (-2*eye(N) + diag(ones(N-1,1),1) + diag(ones(N-1,1),-1));
A = eye(N) - Dt*epsilon*theta_t*delta2;
Ainv = A^(-1);
B = eye(N) + Dt*epsilon*(1-theta_t)*delta2;
C = Dt*epsilon/(Dx^2)*[u_l;zeros(N-2,1);u_r];

% Initialisation
Un = u0(X_A);
U_out = zeros(N,Nt);
X_D_out = zeros(N,Nt,K);
Unplus1 = Un;
X_D = X_A - Dt*Unplus1;
if track_front
  bigXstar = zeros(Nt+1,1);
end % if track_front

for tt = 1:Nt
    for kk = 1:K
      % Calculate departure points
      X_D_old = X_D;
      X_D = X_A - Dt*Unplus1;
      X_D(X_D<x_l) = x_l;
      X_D(X_D>x_r) = x_r;
      for ll = 1:2
        Un_D = interp1(X_A_bc,[u_l;Un;u_r],X_D);
        X_D = X_A - Dt*(theta_x*Unplus1 + (1 - theta_x)*Un_D);
        X_D(X_D<x_l) = x_l;
        X_D(X_D>x_r) = x_r;
      end
      R_D = interp1(X_A_bc,[u_l;B*Un + (1-theta_t)*C;u_r],X_D); % No u_xx on x_l, x_r
      %ENO_pp = interp_ENO(X_A_bc,[u_l;B*Un + (1-theta_t)*C;u_r]);
      %R_D = ppval(ENO_pp, X_D);
      Unplus1 = Ainv*(R_D + theta_t*C);
      X_D_out(:,tt,kk) = X_D;
      X_D_diff = X_D - X_D_old;
    end % for kk
    Un = Unplus1;
    U_out(:,tt) = Unplus1;
    if track_front
      [~,bigXstar(tt+1)] = get_m_x(Un,X_A,c,alpha_0);
    end % if track_front
end % for tt

if plotting
  % Plot u0(X) and U(plottimes) whilst holding
  h1 = figure;
  set(h1,'defaulttextinterpreter','latex');
  x_hr = linspace(x_l,x_r,1000);
  plot(x_hr, u0(x_hr),'k',0,c,'kx');
  hold on
  xlabel '$x$';
  ylabel '$U(x,t)$';
  for tt = plottimes
    plot(X_A_bc,[u_l;U_out(:,tt);u_r],'k')
    ylim(plotlims)
    drawnow
    % Put a cross on the front
    if track_front
      plot(bigXstar(tt+1),c,'kx')
    end
  end

%  myt = 1001;
%  for kk = 1:K
%    for i = 1:N
%      plot([X_D_out(i,myt,kk),X_A(i)],[0,1]), hold on
%    end % for i
%    pause, hold off
%  end % for kk
end % if plotting

path(old_path);
end % function
