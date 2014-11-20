function ut = fwd_euler(uin, xin, epsilon)
% function giving the time derivative of viscous burgers equation
%  u_t = -u*u_x + epsilon * u_(xx),
% using a finite difference spacial discretisation, with dirichlet 
% boundary conditions. Could even be used for implicit RK methods.

% For use in burgersSLMM.m for the forward euler approximation to
% u.

N = length(uin);

ut = zeros(size(uin));
Dx = diff(xin);

ut(1) = 0;
ut(N) = 0;

for ii = 2:(N-1)
    ut(ii) = (uin(ii-1) + uin(ii) + uin(ii+1))/3 * ...
        (uin(ii+1) - uin(ii-1))/(Dx(ii)+Dx(ii-1)) + ...
        epsilon*(...
        (uin(ii+1)-uin(ii))/(1/2*Dx(ii)*(Dx(ii) + Dx(ii-1))) - ...
        (uin(ii)-uin(ii-1))/(1/2*Dx(ii-1)*(Dx(ii) + Dx(ii-1))));
end

