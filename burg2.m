function [Un,X_A] = burg2(N,Nt)
% A light-weight SISL burgers solver.
%
% Uniform grid, linear interpolation, Crank-Nicholson-like Departure points

%clear all

% Parameters
%N = 300;
%Nt = 3000;
K = 40;

t0 = 0;
tmax = 1.5;
epsilon = 0.0001;
alpha_0 = 0.5;

%x_l = 0;
%x_r = 1;
%u0 = @(x) sin(2*pi*x) + 1/2*sin(pi*x) - 0.5;
%plotlims = [-2.5 2.5];
x_l = -1;
x_r = 4;
c = 1;
alpha_ = 0.1;
u0 = @(x) c - alpha_*tanh(alpha_/(2*epsilon)*(x - c*t0));

%tmax = 15;
%x_l = -10;
%x_r = 40;
%u0 = @(x) c - alpha_*tanh(alpha_/(2*epsilon)*(x - c*t0));


plotting = 0;
plotlims = [0 1.6];

% Setup
Dx = (x_r-x_l)/(N+1);
Dt = tmax/Nt;
u_l = u0(x_l);
u_r = u0(x_r);

X_A_bc = (x_l:Dx:x_r)';
X_A = X_A_bc(2:end-1);

delta2 = 1/(Dx^2) * ...
    (-2*eye(N) + diag(ones(N-1,1),1) + diag(ones(N-1,1),-1));
A = eye(N) - Dt*epsilon*alpha_0*delta2;
Ainv = A^(-1);
B = eye(N) + Dt*epsilon*(1-alpha_0)*delta2;
C = Dt*epsilon/(Dx^2)*[u_l;zeros(N-2,1);u_r];

% Initialisation
Un = u0(X_A);
U_out = zeros(N,Nt);
X_D_out = zeros(N,Nt,K);
Unplus1 = Un;
X_D = X_A - Dt*Unplus1;

for tt = 1:Nt
    for kk = 1:K
      % Calculate departure points
      X_D_old = X_D;
      X_D = X_A - Dt*Unplus1;
      X_D(X_D<x_l) = x_l;
      X_D(X_D>x_r) = x_r;
      for ll = 1:2
        Un_D = interp1(X_A_bc,[u_l;Un;u_r],X_D);
        X_D = X_A - Dt*(0.5*Unplus1 + 0.5*Un_D);
        X_D(X_D<x_l) = x_l;
        X_D(X_D>x_r) = x_r;
      end
      R_D = interp1(X_A_bc,[u_l;B*Un + (1-alpha_0)*C;u_r],X_D); % No u_xx on x_l, x_r
      %ENO_pp = interp_ENO(X_A_bc,[u_l;B*Un + (1-alpha_0)*C;u_r]);
      %R_D = ppval(ENO_pp, X_D);
      Unplus1 = Ainv*(R_D + alpha_0*C);
      X_D_out(:,tt,kk) = X_D;
      X_D_diff = X_D - X_D_old;
    end % for kk
    Un = Unplus1;
    U_out(:,tt) = Unplus1;
end % for tt

if plotting
  for tt = 1:Nt
  plot(X_A_bc,[u_l;U_out(:,tt);u_r])
  ylim(plotlims)
  drawnow
  end

%  myt = 1001;
%  for kk = 1:K
%    for i = 1:N
%      plot([X_D_out(i,myt,kk),X_A(i)],[0,1]), hold on
%    end % for i
%    pause, hold off
%  end % for kk
end % if plotting

end % function
