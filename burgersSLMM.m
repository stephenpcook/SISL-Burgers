function [U_A,X_An1,bigXstar,bigDxMin] = burgersSLMM(N, tN, param_file)
%BURGERSSLMM Solving Burgers with a SL method with moving meshes.
%
% [U_A,X_An1] = BURGERSSLMM(N, tN, param_file) returns the end-time
% numerical solution and the arrival mesh for N points and tN time steps,
% with all other parameters defined with the .mat filename given in
% param_file.
%
% [U_A,X_An1] = BURGERSSLMM(N, tN) is the same but only uses the default
% parameter values from options/param_defaults.mat, as defined by
% gen_param_defaults.
%
% [U_A,X_An1,bigXstar] = BURGERSSLMM(N, tN, ...) also makes a call to
% get_m_x and finds the points x_star such that
% U(x_star) = (u(x_l) + u(x_r))/2.
%
% [U_A,X_An1,bigXstar,bigDxMin] = BURGERSSLMM(N, tN, ...) also returns the
% minimum mesh spacing over time.
%
% This uses the following external files:
%  get_m_x.m
%  mmpde5.m
%  fwd_euler.m
% And from subdirectories
%  ~/dos/MATLAB/mk_video.m
%  ./interpolation/interp1cubicL.m
%  ./interpolation/ppval_lim.m
%  ./mm_suite/Eqd1dExact.m
%
% See also: BURG2 GEN_PARAM_DEFAULTS GET_M_X MMPDE5 FWD_EULER INTERP1CUBICL
% PPVAL_LIM EQD1DEXACT

% Author: Stephen P. Cook <s.cook@bath.ac.uk>
% Date: 03-01-2018
old_path = addpath([pwd,'\interpolation'],...
                   [pwd,'\mm_suite'],...
                   [pwd,'\options']);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%           Default Parameters           %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% General Parameters
%%N=21;       % Level of spacial discretisation
%%tN = 150;     % Number of timesteps
%K=10;         % Departure point iterations
%theta_t=1/2;  % Theta for the Theta method in time.
%theta_x=1/2;  % Theta for the Theta-method for departure points.
%epsilon=0.0001 ; % Epsilon in the PDE
%
%%% Domain Parameters
%x_l = -1;
%x_r = 4;
%
%t0 = 0;
%tmax = 1.5;   % Final time
%
%%% Initial and boundary conditions
%%u0 = @(x) (sin(2*pi*x) + 1/2*sin(pi*x)) + 1;
%%u_l = 1; u_r = 1;
%c = 1;
%alpha_0 = 0.1;
%u0 = @(x) c - alpha_0*tanh(alpha_0/(2*epsilon)*(x - c*t0));
%u_l = u0(x_l);
%u_r = u0(x_r);
%
%%% Plot parameters
%%plotlims = [-0.5, 1.5];
%plotlims = [c-alpha_0-0.1, c+alpha_0+0.1];
%plotting = 1;
%
%%% Mesh Parameters
%%mesh_movement = 'static';
%%mesh_movement = 'prescribed';
%mesh_movement = 'moving-exact';  % Mesh movement type
%%mesh_movement = 'moving-relax';  % Mesh movement type
%limiter = 1;        % Flux limiter for interpolation
%interpolation = 'linear';
%%interpolation = 'CSpline';
%%interpolation = 'CLagrange';
%%interpolation = 'ENO';
%%interpolation = 'pchip';
%
%%% Monitor function parameters
%b = 0.1;
%m=@(x,u,u_x, u_xx) sqrt(b + u_x.^2);
%%m=@(x,u,u_x,u_xx)ones(size(u))
%p_smooth = 5;
%tau =1;
%with_euler = 1;
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%            START OF CODE              %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

load('params_default');
if nargin==3
  % There could be an issue here if we try to define an anonymous function
  % in a .mat file which is not in the path. Be aware
  load(param_file);
end

if nargout>=3
  track_front = 1;
else
  track_front = 0;
end
if nargout==4
  track_min_dx = 1;
else
  track_min_dx = 0;
end

plotting = 0;

Dx= (x_r-x_l)./(N+1); % Space step
Dt = (tmax-t0)/tN;  % Time step

% Matrices to be used. Tridiagonal, so could be solved more efficiently

% delta2 - \delta_x^2, a tridiagonal matrix, finite difference second
% derivative operator.
delta2 = (-2*eye(N) + diag(ones(N-1,1),1) + diag(ones(N-1,1),-1))./(Dx^2);
% LHS and RHS of SISL formulation
M_RHS = eye(N) + Dt * (1 - theta_t) * epsilon * delta2;
M_LHS = eye(N) - Dt * theta_t * epsilon * delta2;
BC = [u_l/(Dx^2);zeros(N-2,1);u_r/(Dx^2)];

% Initialisation
X = x_l + (1:N)'*Dx;
Dxi = 1/(N+1);
Un = u0(X);
X_An = X;
TT = t0 + (0:tN)'*Dt;
XX = zeros(length(TT),length(X));
XX(1,:) = X;
dept_convergence = zeros(K,length(TT));
if track_front
  bigXstar = zeros(tN,1);
end % if track_front
bigDxMin = zeros(tN,1);
uout= XX;
uout(1,:) = Un;
X_An1 = X_An;
DX_An = diff([x_l;X_An;x_r]);
DX_An1 = diff([x_l;X_An1;x_r]);

del2n1 = delta2;
del2n = NaN(N);
BCn1 = BC;
BCn = BC;

% Outer loop. Timestep
for tt = 1:length(TT)
    t = TT(tt);
% Timestep initialisation
% Cubic spline of Un (including the constant end points)

% Mesh movement
if strcmp(mesh_movement, 'static')
  vtitle = ['Semi-Lagrangian Burgers, static mesh, N = ',num2str(N)];
else
switch mesh_movement
    case 'prescribed'
        % prescribed mesh, loaded from file, innit.
	% Not yet done, static for now.
        X_An1 = X;
        if tt>=floor(tN/2)
          Xi = (1:N)'./(N+1);
          alpha_p = 0.9;
          X_An1 = x_l + (alpha_p*Xi.^2 + (1-alpha_p)*Xi)*(x_r-x_l);
        end
    case 'moving-exact'
        U = [u_l;Un;u_r];
        X_ = [x_l;X_An;x_r];
        if with_euler
            Uhat = U + Dt*fwd_euler(U,X_,epsilon);
            DUhatDX = [Uhat(2)-Uhat(1);diff(Uhat)]./[X_(2)-X_(1);diff(X_)];
            D2UhatDX2 = [0;del2n1*Uhat(2:end-1) + BCn1;0];
            D2UhatDX2(1) = D2UhatDX2(2);
            D2UhatDX2(end) = D2UhatDX2(end-1);
            M = m(X_,Uhat,DUhatDX,D2UhatDX2);
        else
            DUDX = [Un(1) - u_l;diff(U)]./[X_An(1)-x_l;diff(X_)];
            D2UDX2 = [0;del2n1*Un + BCn1;0];
            D2UDX2(1) = D2UDX2(2);
            D2UDX2(end) = D2UDX2(end-1);
            M = m(X_,U,DUDX,D2UDX2);
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
	vtitle = ['Semi-Lagrangian Burgers, moving mesh (exact), N = '...
	,num2str(N)];
    case 'moving-relax'
        % This is different from the above method, in that
        % we solve a moving mesh PDE to get the movement,
        % and to do the we take M to sit at the midpoints.

        X_ = [x_l;X_An;x_r];
        U = [u_l;Un;u_r];
        DUDX = diff(U)./diff(X_);    % ! Sits at the midpoints
        D2UDX2 = del2n1*Un + BCn1;
        M = m(X_,U,DUDX,D2UDX2); % This might require changing for
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
        M=M./(trapz(1/2*([x_l;X_An]+[X_An;x_r]),M)); % normalise
        M = 1/2 + M./2;             % Average
        %
        x = X_;
        % Spline it!
        M_pp = spline(x,mmpde5(M,x,Dxi,tau));
        [~,ode_out] = ode15s(@(t,x_)ppval(M_pp,x_),[0 Dt],x);
        X_An1 = ode_out(end,:)';
        X_An1 = X_An1(2:end-1);
	vtitle = ['Semi-Lagrangian Burgers, moving mesh (relax), N = '...
	,num2str(N)];
    otherwise
        disp('Unknown movement type, using static')
        X_An1 = X;
end % switch

DX_An = diff([x_l;X_An;x_r]);
DX_An1 = diff([x_l;X_An1;x_r]);

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
rhs0 = u_l;
rhsN = u_r;

rhs_A = [rhs0; (M_RHS * Un + Dt*(1-theta_t)*epsilon*BCn) ; rhsN];

switch interpolation
    case 'linear'
        f_rhs = griddedInterpolant([x_l;X_An;x_r],rhs_A, 'linear');
        f_Un = griddedInterpolant([x_l;X_An;x_r],[u_l;Un;u_r],'linear');
    case 'CSpline'
        % TODO Not sure if this is what we want
        f_rhs = griddedInterpolant([x_l;X_An;x_r],rhs_A, 'spline');
        f_Un = griddedInterpolant([x_l;X_An;x_r],[u_l;Un;u_r],'spline');
    case 'CLagrange'
        f_rhs = griddedInterpolant([x_l;X_An;x_r],rhs_A, 'cubic');
        f_Un = griddedInterpolant([x_l;X_An;x_r],[u_l;Un;u_r],'cubic');
    case 'ENO'
        pp_rhs = interp_ENO([x_l;X_An;x_r], rhs_A);
        pp_f = interp_ENO([x_l;X_An;x_r],[u_l; Un ; u_r]);
        f_rhs = @(Xq) ppval(pp_rhs,Xq);
        f_Un = @(Xq) ppval(pp_f,Xq);
    case 'pchip'
        f_rhs = griddedInterpolant([x_l;X_An;x_r],rhs_A, 'pchip');
        f_Un = griddedInterpolant([x_l;X_An;x_r],[u_l;Un;u_r],'pchip');
    case 'hermite'
        if limiter
            f_rhs = @(Xq)interp_hermite_lim(...
                [x_l;X_An;x_r],rhs_A,Xq,'hyman');
            f_Un = @(Xq)interp_hermite_lim(...
                [x_l;X_An;x_r],[u_l;Un;u_r],Xq,'hyman');
        else
            D_rhs = calc_gradients([x_l;X_An;x_r],rhs_A,'hyman');
            f_rhs = @(Xq)eval_hermite([x_l;X_An;x_r],rhs_A,D_rhs,Xq);

            D_Un = calc_gradients([x_l;X_An;x_r],[u_l;Un;u_r],'hyman');
            f_Un = @(Xq)eval_hermite([x_l;X_An;x_r],[u_l;Un;u_r],D_Un,Xq);
        end
end

% Initial guess of departure points.
X_D = X_An1 - Dt*Un;
X_D(X_D<x_l) = x_l;
X_D(X_D>x_r) = x_r;

% Evaluate RHS at the departure points.
rhs_D = f_rhs(X_D);
if limiter
  rhs_D = feval_lim(X_D,rhs_D,[x_l;X_An;x_r],rhs_A);
end % if limiter

 % Inner loop.
 for k = 1:K

 % Solve the implicit equation (Thomas algorithm implemented later)
 U_A = M_LHS\(rhs_D + Dt*theta_t*epsilon*BCn1);
 % New guess at the departure points (theta-method) and recalculate
 % RHS_D and U_A.
 for ll = 1:2
   Un_D = f_Un(X_D);
   if limiter
     Un_D = feval_lim(X_D,Un_D,[x_l;X_An;x_r],[u_l;Un;u_r]);
   end % if limiter
   X_D_old = X_D;
   X_D = X_An1 - Dt*(theta_x*U_A + (1-theta_x)*Un_D);
   X_D(X_D<x_l) = x_l;
   X_D(X_D>x_r) = x_r;
 end % for ll
 dept_convergence(k,tt) = norm(X_D - X_D_old);

 % Evaluate RHS at the departure points.
 rhs_D = f_rhs(X_D);
 if limiter
   rhs_D = feval_lim(X_D,rhs_D,[x_l;X_An;x_r],rhs_A);
 end % if limiter
 end % for k

 U_A = M_LHS\(rhs_D + Dt*theta_t*epsilon*BCn1);
 % End of Inner Loop


%%% Prep for next timestep. %%%
X_An = X_An1;

Un = U_A;
uout(tt,:) = Un;
if track_front
  [~,bigXstar(tt)]=get_m_x(U_A,X_An1,c, alpha_0);
end % if track_front
if track_min_dx
  bigDxMin(tt) = min(diff([x_l;X_An1;x_r]));
end % if track_min_dx
%%% And plot %%%
if plotting
  plot([x_l;X_An1;x_r],[u_l;U_A;u_r])
  title(['t = ',num2str(t)]), ylim(plotlims)
  drawnow()
  XX(tt,:) = X_An;
  bigX_D(tt,:) = X_D;
end % if plotting
end % for tt

if plotting
% pause
 h1=figure;
 set(h1,'defaulttextinterpreter','latex');
 for i = 1:(floor(N/25)):N
    plot(XX(:,i),TT,'k')
    hold on
 end % for i
 hold off
 %title 'Mesh trajectories (showing 26 mesh points)'
 xlabel('$X_A$','FontSize',18)
 ylabel('$t$','FontSize',18)
 h2=figure;
 set(h2,'defaulttextinterpreter','latex');
 semilogy(TT,min(diff(XX')),'k')
 %title 'Minimum mesh spacing over time'
 xlabel('$t$','FontSize',18)
 ylabel('$\displaystyle{\min_j(X_j(t) - X_{j-1}(t))}$','FontSize',18)
end % if plotting

path(old_path);
end % function main
