function [un,x] = burgersTom(nx,nt)
% SISL implmentation of burgers equation

%nx = 81;
%nt = 150;
%dt = 0.01;
dt = 1.5/nt;
alpha_0 = 0.5;
epsilon = 10^-4;
c = 1;
h = 0.5;
n_dept = 2;
n_outer = 2;

plotting = 0;

x = zeros(nx,1);
for i = 1:nx
  x(i) = 5*(i-1)/(nx-1);
end
dx = x(2)-x(1);
un = zeros(nx,1);

% initial conditions
un(:) = c - h;
un(1) = c + h;
unp1 = un;
un_eul = un;

% build lhs
lhs = zeros(nx,nx);
for i = 1:nx
  lhs(i,i) = 1;
end
for i = 2:nx-1
  lhs(i,i) =  lhs(i,i) - 2*alpha_0*dt*epsilon/dx^2;
  lhs(i,i+1) = alpha_0*dt*epsilon/dx^2;
  lhs(i,i-1) = alpha_0*dt*epsilon/dx^2;
end

for n = 1:nt
% compute rn
  rn = un;
  rn(2:nx-1) = rn(2:nx-1) - (1-alpha_0)*dt*epsilon*(un(3:nx) - 2*un(2:nx-1) + un(1:nx-2))/dx^2;
  for outer = 1:n_outer
% compute departure points
    xd = x;
    for it = 1:n_dept
      un_d = interp1(x,un,xd,'linear');
      %un_pp = interp_ENO(x,un);
      %un_d = ppval(un_pp, xd);
      xd = x - dt/2*(unp1 + un_d);
      for i = 1:nx
        if ( xd(i) < 0 )
	  xd(i) = 0;
	elseif ( xd(i) > 5 )
	  xd(i) = 5;
	end
      end
    end
% interpolate rn
    rn_d = interp1(x,rn,xd,'linear');
    %rn_pp = interp_ENO(x,rn);
    %rn_d = ppval(rn_pp, xd);
% solve
    unp1 = lhs\rn_d;

  end
  un = unp1;
% Eulerian solution
  unp1_eul = un_eul;
  unp1_eul(2:nx-1) = unp1_eul(2:nx-1) - dt*epsilon*(un_eul(3:nx) - 2*un_eul(2:nx-1) + un_eul(1:nx-2))/dx^2;
  unp1_eul(2:nx-1) = unp1_eul(2:nx-1) - dt*1/6*(un_eul(2:nx-1)+4*un_eul(2:nx-1)+un_eul(2:nx-1)).*(un_eul(3:nx) - un_eul(1:nx-2))/(2*dx);
  un_eul = unp1_eul;
end

u_ex = c - h*tanh(h/(2*epsilon)*(x - c*dt*nt));

if plotting
    %figure; plot(x,unp1,'b',x,unp1_eul,'g',x,u_ex,'k','linewidth',2);
    figure; plot(x,unp1,'b',x,u_ex,'k','linewidth',2);
    ylim([c - h - 0.1, c + h + 0.1]);
end % if plotting

end % function
