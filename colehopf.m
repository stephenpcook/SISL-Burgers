% u_t + (u^2/2)_x = epsilon*u_{xx}
% u(0,t) = u(1,t) = 0
% u(x,0) = u_0(x).

% rho(x,t) = sum_{n=1}^{infinity} A_n * X_n(x) * T_n(t)
% X_n(x) = sin(n * pi * x)
% T_n(t) = exp(-epsilon * n^2 * pi^2 * t)
% A_n = 2*int_0^1[exp(-1/(2*epsilon)*int_0^x[u_0(s)]ds)*cos(n*pi*x)]dx

epsilon = 1/200;
N = 150;

u_0 = @(x)sin(2*pi*x);
w_0 = @(x)1/pi*(sin(pi*x)).^2;

A0 = 2 * quad(@(x)exp(-1/(2*epsilon)*w_0(x)),0,1);
for n = 1:N
    A(n) = 2 * quad(@(x)exp(-1/(2*epsilon)*w_0(x)).*cos(n*pi*x),0,1);
end

x = (0:200)'./200;
numerator = zeros(size(x));
denominator = A0/2 * ones(size(x));
for n=1:N
     %rho_0 = rho_0 + A(n)*cos(n*pi*x);
     numerator= numerator + 2*epsilon*pi*A(n)*n*sin(n*pi*x);
     denominator= denominator + A(n)*cos(n*pi*x);
end
u = numerator./denominator;
plot(u)