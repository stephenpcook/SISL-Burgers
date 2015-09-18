% Interpolation error of moving meshes for a prescribed travelling wave.

clf

% Params
epsilon = 0.01;
b = 1;
t_N = 100;
%t_max = 2;
x_N = 100;
x_c_N = 11;
x_l = 1;
x_r = 0;
alpha_0 = (x_l - x_r)/2;

c = (x_l + x_r)/2;
t_max = 1/c;

beta = alpha_0/(2*epsilon);
T = t_max*(0:(t_N-1))./(t_N-1);
dt = t_max/t_N;
dx = 1/x_N;



u = @(x,t) c - alpha_0*tanh(beta*(x-c*t));
u_x = @(x,t) -alpha_0*beta*(1-tanh(beta*(x-c*t)).^2);
u_xx = @(x,t) 2*alpha_0*beta^2*(tanh(beta*(x-c*t)) - ...
                               tanh(beta*(x-c*t)).^3);
x = (0:(x_N-1))./(x_N-1);
x_coarse = (0:(x_c_N-1))./(x_c_N-1);

M.type = 'points';
M.x = x;

xout = zeros(t_N,x_c_N);
x_D = xout;
error = xout;
error_s = xout;

for ii = 1:t_N
    tt = T(ii);
    U = u(x,tt);
    % Mesh movement
    dUdX = u_x(x,tt);
    M.M = sqrt(b + dUdX.^2);
    xout(ii,:) = Eqd1dExact(x_coarse,M);

    subplot(1,2,1)
    plot(x,u(x,tt))
    drawnow()
end

pause

for ii = 1:(t_N-1)
    % departure points
    x_D(ii,:) = xout(ii+1,:) - dt*...
    1/2*(u(xout(ii+1,:),T(ii+1)) + u(xout(ii,:),T(ii)));

end

F = @(x,t) u(x,t) + 1/2*epsilon*dt*u_xx(x,t);
for ii = 1:(t_N-1)
    % interpolation error
    % We are interpolating
    %   u + epsilon*dt*1/2*u_xx
    %LHS_A = u(xout(ii,:),tt) + 1/2*epsilon*dt*u_xx(xout(ii,:),tt);
    %LHS_D = interp1(xout(ii,:),LHS_A,x_D(ii,:));
    LHS_D = interp1(xout(ii,:),F(xout(ii,:),T(ii)),x_D(ii,:));
    error(ii,:) = abs(F(x_D(ii,:),T(ii)) - LHS_D);
end
% Get rid of NaN, from extrapolation
error(isnan(error)) = 0;
emax = max(max(error));
for ii = 1:(t_N-1)
    plot(error(ii,:));
    ylim([0,1.1*emax]);
    drawnow()
end

pause

plot(max(error,[],2),T)
xlabel('Maximum interpolation error')
ylabel('t')

subplot(1,2,2)
hold on
for jj = 1:x_c_N
    plot(xout(:,jj),T','b-')
    plot(x_D(1:(end-1),jj),T(1:end-1)','r-')
end
xlabel('Mesh movement (blue) and departure points (red)')
ylabel('t')
hold off

pause
% Now without Mesh movement

for ii = 1:(t_N-1)
    % departure points
    x_D_c(ii,:) = x_coarse - dt*u(x_coarse,T(ii+1));

end
for ii = 1:(t_N-1)
    % interpolation error
    LHS_D = interp1(x_coarse,F(x_coarse,T(ii)),x_D_c(ii,:));
    error_s(ii,:) = abs(F(x_D_c(ii,:),T(ii)) - LHS_D);
end

% Get rid of NaN, from extrapolation
error_s(isnan(error_s)) = 0;
e_s_max = max(max(error_s));
subplot(1,2,1)
for ii = 1:(t_N-1)
    plot(error_s(ii,:));
    ylim([0,1.1*e_s_max]);
    drawnow()
end

pause

plot(max(error_s,[],2),T)
xlabel('Maximum interpolation error')
ylabel('t')

subplot(1,2,2)
hold off
for jj = 1:x_c_N
    plot([x_coarse(jj),x_coarse(jj)],[T(1) T(end)],'b-') % Sort out this line
    hold on
    plot(x_D_c(1:(end),jj),T(1:end-1)','r-')
end
xlabel('Mesh points (blue) and departure points (red)')
ylabel('t')
hold off
