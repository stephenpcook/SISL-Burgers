function [m,x_star] = get_m_x(U,X,c,alpha_0)
% Calculate the final-time gradient, location and width of a tanh profile
%
% [m,x_star] = get_m_x(U,X,c,alpha_0)

kk = 1;
while U(kk)>c
  kk = kk+1;
end % while
kk = kk-1;
m = (U(kk+1) - U(kk))/(X(kk+1) - X(kk));
%m2 = max(abs(diff(U)./diff(X)));
x_star = ((c - U(kk))/m + X(kk));

end % function get_m_x
