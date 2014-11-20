% A light-weight SISL burgers solver.
%
% Uniform grid, linear interpolation, Crank-Nicholson-like Departure points

clear all
N = 101;
Nt = 1500;
K = 4;

Dx = 1/(N+1);
Dt = 2/Nt;
u0 = @(x) sin(2*pi*x) + 1/2*sin(pi*x);
epsilon = 0.002;
alpha = 0.5;

X_A = (Dx:Dx:(1-Dx))';
X_A_bc = (0:Dx:1)';

delta2 = 1/(Dx^2) * ...
    (-2*eye(N) + diag(ones(N-1,1),1) + diag(ones(N-1,1),-1));
A = eye(N) - Dt*epsilon*alpha*delta2;
Ainv = A^(-1);
B = eye(N) + Dt*epsilon*(1-alpha)*delta2;

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
      for ll = 1:2
        Un_D = interp1(X_A_bc,[0;Un;0],X_D);
        X_D = X_A - Dt*(0.5*Unplus1 + 0.5*Un_D);
      end
      R_D = interp1(X_A_bc,[0;B*Un;0],X_D);
      Unplus1 = Ainv*R_D;
      X_D_out(:,tt,kk) = X_D;
      X_D_diff = X_D - X_D_old;
    end % for kk
    Un = Unplus1;
    U_out(:,tt) = Unplus1;
end % for tt

for tt = 1:Nt
plot(X_A_bc,[0;U_out(:,tt);0])
ylim([-2,2])
drawnow
end

myt = 1001;
for kk = 1:K
  for i = 1:N
    plot([X_D_out(i,myt,kk),X_A(i)],[0,1]), hold on
  end % for i
  pause, hold off
end % for kk
