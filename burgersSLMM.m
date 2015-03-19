%function [U_A,X_An1] = burgersSLMM(N,tN)
% Solving Burgers with a SL method with moving meshes.
%
% This uses the following external files.
% ./mmpde5.m
% ./fwd_euler.m
% ~/dos/MATLAB/Eqd1dExact.m
% ~/dos/MATLAB/interp1cubicL.m
% ~/dos/MATLAB/ppval_lim.m
% ~/dos/MATLAB/mk_video.m

clear all
clf

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%               Parameters              %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% General Parameters
N=21;       % Level of spacial discretisation
tN = 150;     % Number of timesteps
K=40;         % Departure point iterations
theta_t=1/2;  % Theta for the Theta method in time.
theta_x=1/2;  % Theta for the Theta-method for departure points.
epsilon=0.0001 ; % Epsilon in the PDE

% Domain Parameters
x0 = -1;
x1 = 4;
Dx= (x1-x0)./(N+1); % Space step

t0 = 0;
tmax = 1.5;   % Final time
Dt = (tmax-t0)/tN;

% Initial and boundary conditions
%u0 = @(x) (sin(2*pi*x) + 1/2*sin(pi*x)) + 1;
%u_l = 1; u_r = 1;
c = 1;
alpha = 0.5;
u0 = @(x) c - alpha*tanh(alpha/(2*epsilon)*(x - c*t0));
u_l = u0(x0);
u_r = u0(x1);

% Plot parameters
%plotlims = [-0.5, 1.5];
plotlims = [c-alpha-0.1, c+alpha+0.1];
plotting = 1;

% Mesh Parameters
%mesh = 'static';
%mesh = 'prescribed';
mesh = 'moving-exact';  % Mesh movement type
%mesh = 'moving-relax';  % Mesh movement type
limiter = 1;        % Flux limiter for interpolation
interpolation = 'linear';
%interpolation = 'CSpline';
%interpolation = 'CLagrange';
%interpolation = 'pchip';

% Monitor function parameters
b = 0.1;
m=@(x,u,uprime) sqrt(b + uprime.^2);
%m=@(x,u,uprime)ones(size(u))
p_smooth = 5;
tau =1;
with_euler = 1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%            START OF CODE              %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Matrices to be used. Tridiagonal, so could be solved more efficiently

% del2 - \delta_x^2, a tridiagonal matrix, finite difference second 
% derivative operator.
del2 = (-2*eye(N) + diag(ones(N-1,1),1) + diag(ones(N-1,1),-1))./(Dx^2);
% LHS and RHS of SISL formulation
M_RHS = eye(N) + Dt/(Dx^2) * (1 - theta_t) * epsilon * del2;
M_LHS = eye(N) - Dt/(Dx^2) * theta_t * epsilon * del2;
BC = [u_l/(Dx^2);zeros(N-2,1);u_r/(Dx^2)];

% Initialisation
X = x0 + (1:N)'*Dx;
Dxi = 1/(N+1);
Un = u0(X);
X_An = X;
TT = t0 + (0:tN)'*Dt;
XX = zeros(length(TT),length(X));
XX(1,:) = X;
dept_convergence = zeros(K,length(TT));
bigXstar = zeros(tN,1);
uout= XX;
uout(1,:) = Un;
%jj=1; % Counting variable for saving entries to uout.
X_An1 = X_An;
DX_An = diff([x0;X_An;x1]);
DX_An1 = diff([x0;X_An1;x1]);

del2n1 = del2;
del2n = eye(N);
BCn1 = BC;
BCn = zeros(N,1);

% Outer loop. Timestep
for tt = 1:length(TT)
    t = TT(tt);
% Timestep initialisation
% Cubic spline of Un (including the constant end points)  

% Mesh movement
if strcmp(mesh, 'static')
else
switch mesh
    case 'static'
        X_An1 = X;
	vtitle = ['Semi-Lagrangian Burgers, static mesh, N = '...
	,num2str(N)];
    case 'prescribed'
        % prescribed mesh, loaded from file, innit.
	% Not yet done, static for now.
        X_An1 = X;
        if tt>=floor(tN/2)
          Xi = (1:N)'./(N+1);
          alpha_p = 0.9;
          X_An1 = x0 + (alpha_p*Xi.^2 + (1-alpha_p)*Xi)*(x1-x0);
        end
    case 'moving-exact'
        %XN = [x0;X_An;x1];
        U = [u_l;Un;u_r];
        X_ = [x0;X_An;x1];
        % TODO what is going on here?
        DUDX = [U(2) - u_l;diff(U)]./[X_An(1)-x0;diff(X_)];
        if with_euler
            Uhat = U + Dt*fwd_euler(U,X_,epsilon);
            DUhatDX = [u_l;diff(Uhat)]./[X_An(1)-x0;diff(X_)];
            M = m(X_,Uhat,DUhatDX);
        else
            M = m(X_,U,DUDX);
        end
        % Smooth
        for iii = 1:p_smooth
            M = [2/3*M(2)+1/3*M(1); M; 2/3*M(end-1) + 1/3*M(end)];
            M = conv(M,[1/4,1/2,1/4],'same');  % Smoothing
            M = M(2:end-1);
        end % iii
        M=M./(trapz(X_,M)); % normalise
        M = 1/2 + M./2;             % Average
        %
        monitor.type = 'points';
        monitor.x = X_;
        monitor.M = M;
        X_An1 = Eqd1dExact(X_, monitor);
        X_An1 = X_An1(2:end-1);
 %       X_An1 = (X_An1 + X_An)./2;
	vtitle = ['Semi-Lagrangian Burgers, moving mesh (exact), N = '...
	,num2str(N)];
    case 'moving-relax'
        % This is different from the above method, in that 
        % we solve a moving mesh PDE to get the movement, 
        % and to do the we take M to sit at the midpoints.
        
        X_ = [x0;X_An;x1];
        %XN = X_;
        U = [u_l;Un;u_r];
        DUDX = diff(U)./diff(X_);    % ! Sits at the midpoints
        M = m(X_,U,DUDX); % This might require changing for
                                  % different monitor functions.
        % Smooth
        for iii = 1:p_smooth
            M = [2/3*M(2)+1/3*M(1); M; 2/3*M(end-1) + 1/3*M(end)];
            % This line is needed to correct the smoothing stencil
            % at the endpoints ( change from [1/4 1/2 1/4] which will
            % use the zero, to [1/3 2/3 0] ).
            M = conv(M,[1/4,1/2,1/4],'same');  % Using convolution
            M = M(2:end-1); % Undo the first line of this loop.
        end % iii
        M=M./(trapz(1/2*([x0;X_An]+[X_An;x1]),M)); % normalise
        M = 1/2 + M./2;             % Average
        %
        x = X_;
        %monitor.M = m([x0;X_An;x1],U,DUDX);
        % Spline it!
        M_pp = spline(x,mmpde5(M,x,Dxi,tau));
        [~,ode_out] = ode15s(@(t,x_)ppval(M_pp,x_),[0 Dt],x);
        X_An1 = ode_out(end,:)';
        %X_An1 = X_ + Dt*mmpde5(M,x,Dxi,tau);
        X_An1 = X_An1(2:end-1);
 %       X_An1 = (X_An1 + X_An)./2;
	vtitle = ['Semi-Lagrangian Burgers, moving mesh (relax), N = '...
	,num2str(N)];
    otherwise
        disp('Unknown movement type, using static')
        X_An1 = X;
end % switch

DX_An = diff([x0;X_An;x1]);
DX_An1 = diff([x0;X_An1;x1]);

del2n = del2n1;
BCn = BCn1;

% Build up the LHS and RHS tridiagonal matrix. Can recycle these 
% and get 2 uses out of each del2 matrix (only redo one)
del2n1(1,1) = -2/(DX_An1(2)*DX_An1(1));
del2n1(1,2) = 2/(DX_An1(2)*(DX_An1(2)+DX_An1(1)));
for ii=2:(N-1)
    del2n1(ii,ii-1) = 2/(DX_An1(ii)*(DX_An1(ii+1)+DX_An1(ii)));
    del2n1(ii,ii) = -2/(DX_An1(ii)*DX_An1(ii+1));
    del2n1(ii,ii+1) = 2/(DX_An1(ii+1)*(DX_An1(ii+1)+DX_An1(ii)));
end
del2n1(N,N-1) = 2/(DX_An1(N)*(DX_An1(N+1)+DX_An1(N)));
del2n1(N,N) = -2/(DX_An1(N+1)*DX_An1(N));

BCn1(1) = 2*u_l/(DX_An1(1)*(DX_An1(1)+DX_An1(2)));
BCn1(N) = 2*u_r/(DX_An1(N+1)*(DX_An1(N)+DX_An1(N+1)));

M_LHS = eye(N) - Dt * theta_t * epsilon * del2n1;

M_RHS = eye(N) + Dt * (1 - theta_t) * epsilon * del2n;
end % if mesh not static

% Make the spline for interpolation
% TODO Could use a higher order term by using u_l, u(1),u(2) and u(3), like in burgers.m
rhs0 = u_l + Dt * (1-theta_t) * epsilon *...
    (DX_An(2)*u_l - (DX_An(1)+DX_An(2))*Un(1) + DX_An(1)*Un(2))/...
    (1/2*DX_An(1)*DX_An(2)*(DX_An(1)+DX_An(2))); % This is an O(h) approx'n to u''
rhsN = u_r + Dt * (1-theta_t) * epsilon *...
    (DX_An(N)*u_r - (DX_An(N)+DX_An(N+1))*Un(N) + DX_An(N+1)*Un(N-1))/...
    (1/2*DX_An(N)*DX_An(N+1)*(DX_An(N)+DX_An(N+1)));
% TODO This is wrong for non-static meshes, this is a slight improvement
% but not much.
rhs0 = u_l;
rhsN = u_r;

switch interpolation
    case 'linear'
        pp_rhs = interp1([x0;X_An;x1],...
          [rhs0; (M_RHS * Un + Dt*(1-theta_t)*epsilon*BCn) ; rhsN],...
            'linear','pp');
    case 'CSpline'
        pp_rhs = spline([x0;X_An;x1],...
        [rhs0; (M_RHS * Un + Dt*(1-theta_t)*epsilon*BCn) ; rhsN]);
    case 'CLagrange'
        pp_rhs = interp1cubicL([x0;X_An;x1],...
          [rhs0; (M_RHS * Un + Dt*(1-theta_t)*epsilon*BCn) ; rhsN]);
    case 'pchip'
        pp_rhs = pchip([x0;X_An;x1],...
          [rhs0; (M_RHS * Un + Dt*(1-theta_t)*epsilon*BCn) ; rhsN]);
end

% Initial guess of departure points.
X_D = X_An1 - Dt*Un;
 % Inner loop. 
 for k = 1:K
 % Evaluate RHS at the departure points.
 if limiter
     rhs_D = ppval_lim(pp_rhs,X_D);
 else
     rhs_D = ppval(pp_rhs, X_D);
 end
 
 % Solve the implicit equation (Thomas algorithm implemented later)
 U_A = M_LHS\(rhs_D + Dt*theta_t*epsilon*BCn1);
 % New guess at the departure points (theta-method) and recalculate 
 % RHS_D and U_A. 
 X_D_old = X_D;
 X_D = X_An1 - Dt*(theta_x*U_A + (1-theta_x)*Un);
 dept_convergence(k,tt) = norm(X_D - X_D_old);
 %rhs_D = ppval(pp_rhs, X_D);
 if limiter
     rhs_D = ppval_lim(pp_rhs,X_D);
 else
     rhs_D = ppval(pp_rhs, X_D);
 end % if limiter
 end % for k
 U_A = M_LHS\(rhs_D + Dt*theta_t*epsilon*BCn1);
 % End of Inner Loop


%%% Prep for next timestep. %%%
X_An = X_An1;

Un = U_A;
%jj=jj+1;
uout(tt,:) = Un;
%[~,bigXstar(tt)]=get_m_x(U_A,X_An1,c);
%%% And plot %%% 
if plotting
  plot([x0;X_An1;x1],[u_l;U_A;u_r])
  title(['t = ',num2str(t)]), ylim(plotlims)
  %if isequal(mesh,'moving')
  %    hold on
  %    plot(monitor.x,monitor.M,'g-')
  %    hold off
  %end
  drawnow()
  XX(tt,:) = X_An;
  bigX_D(tt,:) = X_D;
end % if plotting
end % for t
%pause
figure
for i = 1:N
   plot(XX(:,i),TT)
   hold on
   %plot(bigX_D(:,i),TT,'g-')
end

% Okay, only want ~200 frames in the movie.
t_skip = ceil(length(TT)/200);
mk_video('test.avi',vtitle,TT(1:t_skip:end),uout(1:t_skip:end,:),...
    XX(1:t_skip:end,:))

%end % function main
