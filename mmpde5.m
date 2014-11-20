function xt = mmpde5(Min, xin, Dxi, tau)
% Function for finding the time derivative of MMPDE5
%  x_t = 1/tau * (x_xi * M)_xi,
% with a finite difference approximation in space.
% M should be one point shorter than x, since it sits on half-
% gridpoints.

N = length(xin);
xt = zeros(size(xin));

for ii=2:(N-1)
    xt(ii) = 1/tau*((xin(ii+1) - xin(ii))*Min(ii) - ...
        (xin(ii) - xin(ii-1))*Min(ii-1))...
        /(Dxi^2);
end
